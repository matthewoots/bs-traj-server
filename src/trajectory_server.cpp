/*
* trajectory_server.cpp
*
* ---------------------------------------------------------------------
* Copyright (C) 2022 Matthew (matthewoots at gmail.com)
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License
*  as published by the Free Software Foundation; either version 2
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
* ---------------------------------------------------------------------
*/

#include "trajectory_server.h"

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

namespace trajectory_server
{
    /** @brief Initialize the parameters for the Bspline server
     * without starting the timer
    **/
    void bspline_server::init_bspline_server(
        int _order, double _duration_secs, double _command_interval, int _knot_div) 
    {
        tu.test_chronos();
        command_interval = _command_interval;
        knot_div = _knot_div;
        order = _order;

        // Re-adjust duration so that our knots and divisions are matching
        duration_secs = bsu.get_corrected_duration(
            command_interval, _duration_secs);

        knot_size = bsu.get_knots_size(
            command_interval, duration_secs, knot_div);

        knot_interval = duration_secs / knot_size;

        acceptable_cp_size = knot_size + order;
        std::cout << "[tserver]" <<
            " acceptable control_point " << KGRN << 
            acceptable_cp_size << KNRM << " knot_size " << KGRN <<
            knot_size << KNRM << " knot_interval " << KGRN <<
            knot_interval << KNRM << std::endl;
    }

    /** @brief Start the Bspline timer and saves the start time **/
    void bspline_server::start_bspline_time() 
    {
        stime = std::chrono::system_clock::now();

        // convert duration from secs to milliseconds
        long int duration_msecs = (long int)(duration_secs * 1000);
        // initialize the very first pair of time span (start and end)
        time_span_chronos.clear();
        time_span_chronos.push_back(stime);
        time_span_chronos.push_back(
            stime + milliseconds(duration_msecs)); // from the alias std::chrono::seconds
        
        timespan.clear();
        timespan.push_back(0.0);
        timespan.push_back(duration<double>(time_span_chronos[1] - time_span_chronos[0]).count());
    
        original_timespan.clear();
        original_timespan.push_back(0.0);
        original_timespan.push_back(duration<double>(time_span_chronos[1] - time_span_chronos[0]).count());
    }

    /** @brief Get the corrected control points that fits the knot span
     * @param cp raw control point data that has not yet been cropped to the
     * acceptable size
    **/
    vector<Eigen::Vector3d> bspline_server::get_valid_cp_vector(vector<Eigen::Vector3d> cp)
    {        
        std::cout << "[tserver]" <<
            " cp_size/acceptable " <<
            cp.size() << "/" << acceptable_cp_size << std::endl;
        
        // Crop the control points to what can be given to the bspline
        vector<Eigen::Vector3d> acceptable_cp;
        // If the control point size is bigger, crop to make it smaller
        if (cp.size() > acceptable_cp_size)
            for (int i = 0; i < acceptable_cp_size; i++)
                acceptable_cp.push_back(cp[i]);

        // If the control point size is smaller, we have to change the duration_secs
        else
        {
            acceptable_cp = cp;
            duration_secs = (cp.size() - order) * knot_interval;
        }
            
        return acceptable_cp;
    }
    

    bspline_server::pva_cmd bspline_server::update_get_command_by_time()
    {
        std::lock_guard<std::mutex> path_lock(bs_path_mutex);
        std::lock_guard<std::mutex> cmd_lock(cmd_mutex);
        std::lock_guard<std::mutex> time_lock(time_mutex);

        // Bspline is updated before we reach here;
        // bs_control_points = cp;

        time_point<std::chrono::system_clock> now_time = 
            system_clock::now();
        bspline_server::pva_cmd pva;

        double rel_now_time = duration<double>(now_time - stime).count();
        if ((timespan[1] - rel_now_time) < 0)
        {
            pva.t = -1.0;
            return pva;
        }

        bspline_trajectory::bs_pva_state_3d pva3;
        pva3 = bsu.get_single_bspline_3d(
            order, timespan, bs_control_points, rel_now_time);

        pva.p = pva3.pos[0];
        pva.t = rel_now_time;
        
        // if (!pva3.vel.empty())
        pva.v = pva3.vel[0];
        double _norm = sqrt(pow(pva.v.x(),2) + pow(pva.v.y(),2));
        double _norm_x = pva.v.x() / _norm;
        double _norm_y = pva.v.y() / _norm;

        if (!pva3.acc.empty())
            pva.a = pva3.acc[0];
        
        pva.yaw = atan2(_norm_y,_norm_x);

        return pva;
    }

    /** @brief Get the uniform distributed control point vector from 
     * RRT points and the maximum velocity **/
    vector<Eigen::Vector3d> bspline_server::get_redistributed_cp_vector(
        Eigen::Vector3d current_target_cp,
        vector<Eigen::Vector3d> cp, double max_vel)
    {
        vector<double> temp_time_span;
        if (original_timespan.empty())
        {
            temp_time_span.push_back(0.0);
            temp_time_span.push_back(duration_secs);
        }
        else
            temp_time_span = original_timespan;
        
        double knot_span = bsu.get_dt(
            acceptable_cp_size, order, temp_time_span);
        std::cout << "[bspline_server] get_redistributed_cp_vector knot_span: " <<
            KBLU << knot_span << KNRM << std::endl;
        
        return ctt.uniform_distribution_of_cp(
            current_target_cp, cp, max_vel, knot_span);
    }
    
    /** @brief After optimization or redistribution, update with the Bspline 
     * control points that are found, and update the next timespan that will be used 
    **/
    void bspline_server::update_timespan_and_control_points(
        vector<Eigen::Vector3d> control_points)
    {
        std::lock_guard<std::mutex> time_lock(time_mutex);

        vector<double> previous_timespan;
        if (!timespan.empty()) 
            previous_timespan = timespan;
        else
        {
            previous_timespan.push_back(0.0);
            previous_timespan.push_back(duration_secs);

            std::lock_guard<std::mutex> path_lock(bs_path_mutex);
        
            bs_control_points.clear();
            bs_control_points = control_points;
            for (int i = 0; i < bs_control_points.size(); i++)
                std::cout << "[main_server] bs_control_points : " << KRED << 
                    bs_control_points[i].transpose() << KNRM << std::endl;
            return;
        }

        timespan.clear();
        double end_time_offset = 0.0;
        if (duration_secs < previous_timespan[1]-previous_timespan[0])
            end_time_offset = 
                (previous_timespan[1]-previous_timespan[0]) - duration_secs;


        // [MATLAB] overlap = 2;
        // [MATLAB] cp_used = order + 2;
        // [MATLAB] timespan2 = [timespan1(1) + (cp_used-order) * dt, 
        // timespan1(2) + (cp_used-order) * dt];
        timespan.push_back(
            previous_timespan[0] + (current_cp_idx-order) * knot_interval);
        timespan.push_back(
            previous_timespan[1] - end_time_offset 
            + (current_cp_idx-order) * knot_interval);
    
        std::lock_guard<std::mutex> path_lock(bs_path_mutex);
        
        bs_control_points.clear();
        bs_control_points = control_points;
    }

    /** @brief Return the query control point and update the overlap control 
     * point vectors **/
    Eigen::Vector3d bspline_server::get_current_cp_and_overlap(
        double additional_secs)
    {
        // [MATLAB] overlap = 0;
        // [MATLAB] cp_used = order + 2;
        // [MATLAB] timespan2 = [timespan1(1) + (cp_used-order) * dt, timespan1(2) + (cp_used-order) * dt];
        // [MATLAB] % overlap_size = length(ctrlpt1) - (overlap-order);
        // [MATLAB] for x=1:axis
        // [MATLAB]     r = rand(1,(points+order)-order-overlap);
        // [MATLAB]     ctrlpt2(x,:) = numrange * r; % multiply the random [0-1] by the number range
        // [MATLAB] end

        // [MATLAB] ctrlpt2 = [ctrlpt1(:,cp_used-order+1:cp_used+overlap) ctrlpt2];
    
        double running_time = get_running_time();
        std::cout << "[bspline_server] running_time for current_cp: " 
            << KBLU << running_time << KNRM << "s" << std::endl;

        int current_cp_idx = -1;
        Eigen::Vector3d nan_vector = Eigen::Vector3d(
            std::numeric_limits<double>::quiet_NaN(), 
            std::numeric_limits<double>::quiet_NaN(), 
            std::numeric_limits<double>::quiet_NaN());
        
        if (timespan.empty())
            return nan_vector;
            
        for (int i = 0; i <= knot_size; i++)
        {
            double time_at_knot = timespan[0] + i * knot_interval;
            
            if (running_time - time_at_knot <= 0)
            {
                current_cp_idx = i;
                std::cout << "[bspline_server] idx: " <<
                    KBLU << current_cp_idx << KNRM << " time_at_knot: " <<
                    KBLU << time_at_knot << KNRM << "s" << std::endl;
                break;
            }
        }

        if (current_cp_idx - order + 0 < 0)
            return nan_vector;

        if (current_cp_idx < 0)
            return nan_vector;

        int overlap = (int)ceil(additional_secs / knot_interval);
        
        int size_of_overlapping_cp = (current_cp_idx + overlap) - (current_cp_idx-order+1); 

        overlapping_control_points.clear();
        for (int i = 0; i < size_of_overlapping_cp; i++)
        {
            overlapping_control_points.push_back(
                bs_control_points[current_cp_idx - order + i]);
            std::cout << "[main_server] " << KCYN <<
                current_cp_idx - order + i << KNRM << " bs_control_points : " << KCYN << 
                bs_control_points[current_cp_idx - order + i].transpose() << KNRM << std::endl;
        }

        return bs_control_points[current_cp_idx];
    }

    /** @brief Outdated code **/
    /** @brief For getting index rather than getting by time
     * update_get_command_on_path_by_idx and update_bs_path
     * These 3 functions comes together, this is because update_bs_path have to generate
     * the bspline path and the command functions has to search for the idx
    */
    bspline_trajectory::bs_pva_state_3d bspline_server::get_bs_path(
        vector<Eigen::Vector3d> cp)
    {
        bspline_trajectory::bs_pva_state_3d pva3;
        vector<Eigen::Vector3d> acceptable_cp = 
            get_valid_cp_vector(cp);

        bs_control_points = acceptable_cp;

        pva3 = bsu.get_uni_bspline_3d(
            order, timespan, acceptable_cp, knot_div);

        std::cout << "[tserver]" <<
            " path_size: " <<
            pva3.pos.size() << std::endl;

        return pva3;
    }

    /** @brief Outdated code **/
    /** @brief For getting index rather than getting by time
     * update_get_command_on_path_by_idx and update_bs_path
     * These 3 functions comes together, this is because update_bs_path have to generate
     * the bspline path and the command functions has to search for the idx
    */
    void bspline_server::update_bs_path(vector<Eigen::Vector3d> cp)
    {
        std::lock_guard<std::mutex> path_lock(bs_path_mutex);

        pva_state.rts.clear();
        pva_state.pos.clear(); pva_state.vel.clear(); pva_state.acc.clear();
        pva_state = bspline_server::get_bs_path(cp);
    
    }

    /** @brief Outdated code **/
    /** @brief For getting index rather than getting by time
     * update_get_command_on_path_by_idx and update_bs_path
     * These 3 functions comes together, this is because update_bs_path have to generate
     * the bspline path and the command functions has to search for the idx
    */
    bspline_server::pva_cmd bspline_server::update_get_command_on_path_by_idx()
    {
        std::lock_guard<std::mutex> path_lock(bs_path_mutex);
        std::lock_guard<std::mutex> cmd_lock(cmd_mutex);

        time_point<std::chrono::system_clock> now_time = 
            system_clock::now();
        bspline_server::pva_cmd pva;

        bool early_break = false;
        int idx;
        double rel_now_time = duration<double>(now_time - stime).count();
        // Check the time now and compare to knots
        for(int i = 0; i < (int)pva_state.rts.size(); i++)
        {
            // std::cout << pva_state.rts[i] - rel_now_time << std::endl;
            if((pva_state.rts[i] - rel_now_time) > 0)
            {
                idx = i;
                early_break = true;
                break;
            }
        }
        // if no early break, return update t as -1 to reject the command 
        if (!early_break)
        {
            pva.t = -1.0;
            return pva;
        }
        else
        {
            pva.p = pva_state.pos[idx];
            pva.t = pva_state.rts[idx];
            
            // Update the velocity here but we can also calculate yaw here
            // We can do it 2 ways, using velocity or position difference vector
            double _norm_x, _norm_y;
            if (!pva_state.vel.empty())
            {
                pva.v = pva_state.vel[idx];
                double _norm = sqrt(pow(pva.v.x(),2) + pow(pva.v.y(),2));
                _norm_x = pva.v.x() / _norm;
                _norm_y = pva.v.y() / _norm;
            }
            else
            {
                Vector3d _pos_diff = pva_state.pos[idx] - pva_state.pos[idx-1];
                double _norm = sqrt(pow(_pos_diff.x(),2) + pow(_pos_diff.y(),2));
                _norm_x = _pos_diff.x() / _norm;
                _norm_y = _pos_diff.y() / _norm;
            }
            if (!pva_state.acc.empty())
                pva.a = pva_state.acc[idx];
            
            pva.yaw = atan2(_norm_y,_norm_x);

        }

        return pva;
    }
}