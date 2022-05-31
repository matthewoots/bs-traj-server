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
    vector<Eigen::Vector3d> bspline_server::get_valid_cp_vector(vector<Eigen::Vector3d> cp)
    {        
        std::cout << "[tserver]" <<
            " cp_size/acceptable " <<
            cp.size() << "/" << acceptable_cp_size << std::endl;
        
        // Crop the control points to what can be given to the bspline
        vector<Eigen::Vector3d> acceptable_cp;
        for (int i = 0; i < acceptable_cp_size; i++)
            acceptable_cp.push_back(cp[i]);
            
        return acceptable_cp;
    }

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

    void bspline_server::update_bs_path(vector<Eigen::Vector3d> cp)
    {
        std::lock_guard<std::mutex> path_lock(bs_path_mutex);

        pva_state.rts.clear();
        pva_state.pos.clear(); pva_state.vel.clear(); pva_state.acc.clear();
        pva_state = bspline_server::get_bs_path(cp);
    
    }
    
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

    bspline_server::pva_cmd bspline_server::update_get_command_by_time(
        vector<Eigen::Vector3d> cp)
    {
        std::lock_guard<std::mutex> path_lock(bs_path_mutex);
        std::lock_guard<std::mutex> cmd_lock(cmd_mutex);

        bs_control_points = cp;

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
            order, timespan, cp, rel_now_time);

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

    vector<Eigen::Vector3d> bspline_server::get_redistributed_cp_vector(
        vector<Eigen::Vector3d> cp, double max_vel)
    {
        double knot_span = bsu.get_dt(
            acceptable_cp_size, order, original_timespan);
        
        Eigen::Vector3d current_target_cp;


        return ctt.uniform_distribution_of_cp(
            current_target_cp, cp, max_vel, knot_span);
    }

    bool bspline_server::valid_cp_count_check(size_t cp_size)
    {
        return ((int)cp_size > (knot_size + order));
    }

    double bspline_server::get_running_time()
    {
        time_point<std::chrono::system_clock> now = system_clock::now();
        return duration<double>(now - stime).count();
    }
    
}