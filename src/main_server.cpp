/*
* main_server.cpp
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

#include "main_server.h"

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
    /** @brief Use this function wisely, since check and update may cause an infinite loop
    * If a bad data is given to the rrt node and it cannot complete the validity check
    */
    void main_server::check_and_update_search(
        pcl::PointCloud<pcl::PointXYZ>::Ptr obs,  Eigen::Vector3d first_cp)
    {
        std::lock_guard<std::mutex> pose_lock(pose_update_mutex);
        std::lock_guard<std::mutex> search_points_lock(search_points_mutex);

        if (previous_search_points.empty())
        {
            std::cout << KRED <<
                "[main_server] previous_search_points.size() is empty" <<
                KNRM << std::endl;
            return;
        }

        // Check to see whether the new control point and the previous inputs
        // have any pointclouds lying inside
        if (!ru.check_line_validity_with_pcl(
            previous_search_points[0], first_cp, obs_threshold*1.1, obs))
        {
            previous_search_points.clear();
            std::cout << KRED <<
                "[main_server] connecting point is not obs free" <<
                KNRM << std::endl;
            return;
        }

        int last_safe_idx = -1;
        int initial_size = (int)previous_search_points.size();
        for (int i = 1; i < initial_size; i++)
        {
            if (!ru.check_line_validity_with_pcl(
                previous_search_points[i], previous_search_points[i-1], obs_threshold*1.1, obs))
            {
                last_safe_idx = i-1;
                break;
            }
        }

        if (last_safe_idx >= 0)
        {
            for (int i = last_safe_idx + 1; i < initial_size; i++)
                previous_search_points.erase(previous_search_points.end());
        }

        // std::cout << "[main_server] previous_search_points.size(): " << 
        //     KGRN << previous_search_points.size() << KNRM << std::endl;
    }

    /** @brief Construct the search path from RRT search and from its shortened path */
    void main_server::generate_search_path(
        pcl::PointCloud<pcl::PointXYZ>::Ptr obs)
    {
        std::lock_guard<std::mutex> pose_lock(pose_update_mutex);

        double radii;
        radii = search_radius;
        // std::cout << "[main_server] generate_search_path from " <<
        //     KBLU << current_control_point.transpose() << KNRM << 
        //     " to " << 
        //     KBLU << end.transpose() << KNRM << std::endl;

        // Find a RRT path that is quick and stretches to the end point
        vector<Eigen::Vector3d> path = fe_rrt_server.find_rrt_path(
            previous_search_points, obs, current_control_point, end, radii);
        // If path gives an invalid value, execute some form of emergency
        if (path.empty())
        {
            return;
        }

        // Save global_path
        global_search_path.clear();
        global_search_path = path;
    }

    /** @brief Get the local control points from the RRT module (via distribution) */
    vector<Eigen::Vector3d> main_server::get_local_control_points() 
    {
        // Save new global_search_path into previous_search_points
        previous_search_points.clear();
        for (int i = 0; i < (int)global_search_path.size(); i++)
        {
            // Do not add the start point
            if (i == 0)
                continue;

            previous_search_points.push_back(global_search_path[i]);
        }  
        // std::cout << "[main_server] previous_search_points.size(): " <<
        //     KBLU << previous_search_points.size() << KNRM << std::endl;

        return ts.get_redistributed_cp_vector(
            current_control_point,
            previous_search_points, max_vel);
    }
    
    /** @brief Run the whole algorithm to acquire the control points */
    void main_server::complete_path_generation()
    {
        // Start with getting the current_control_point, and checking the overlap 
        current_control_point = 
            ts.get_current_cp_and_overlap();

        if (isnan(current_control_point[0]) || 
            isnan(current_control_point[1]) || 
            isnan(current_control_point[2]))
        {
            // Use the start point if there is error in getting the control point
            current_control_point = start;
        }

        std::cout << "[main_server] starting from [" <<
            KBLU << current_control_point.transpose() << KNRM << "] to ["
            KBLU << end.transpose() << KNRM << "]" << std::endl;

        // Setup the RRT server parameters 
        fe_rrt_server.reset_parameters(
            vector<Eigen::Vector4d>(),
            min_height, max_height, obs_threshold,
            sub_runtime_error, runtime_error);
        std::cout << "[main_server]" << KBLU <<  
            " 1. reset_parameters complete" << KNRM << std::endl;

        std::lock_guard<std::mutex> cloud_lock(cloud_mutex);

        // Check for any collisions with the local map, and update the valid previous search nodes 
        check_and_update_search(local_cloud, current_control_point);
        std::cout << "[main_server]" << KBLU << 
            " 2. check_and_update_search complete" << KNRM << std::endl;

        bool bypass = false;
        time_point<std::chrono::system_clock> 
            global_start_time = system_clock::now();

        // Check to see whether the previous data extents to the end, keep this data for checks and to recycle
        if (!previous_search_points.empty())
        {
            /** @brief With valid search points and there is no collision 
            * from check_and_update_search(), hence the waypoints can be reused */ 
            if ((previous_search_points[previous_search_points.size()-1]
                - end).norm() < 0.0001)
            {
                std::cout << "[main_server]" << KGRN <<
                    " 3. bypass extract direct goal" <<
                    KNRM << std::endl;

                bypass = true;
            }
        }
        
        /** @brief When there is collision or when the search does not extent to the end */ 
        if (!bypass)
        {
            generate_search_path(local_cloud);    
            std::cout << "[main_server]" << KBLU << 
                " 3. generate_search_path complete" << KNRM << std::endl;
            
            // To add tracking of time to clear the previous search points, so that the check_and_update_search
            // function will be stick with the old search points to check for validity 
            previous_search_points_time.clear();  
            double time_offset = 0.0;
            if (ts.get_duration_from_start_time() < pow(10,6))
                time_offset = ts.get_start_time_for_path();
            double total_time = 0.0 + time_offset;     
            double distance_span = max_vel * ts.get_knot_interval();        
            for (int i = 1; i < global_search_path.size(); i++)
            {
                total_time += 
                    (ceil((global_search_path[i] - global_search_path[i-1]).norm() / distance_span)) * ts.get_knot_interval();
                previous_search_points_time.push_back(total_time);
            }

            // Find and update with newly found distributed control points
            distributed_control_points.clear();
            distributed_control_points = get_local_control_points();
            std::cout << "[main_server]" << KBLU << 
                " 4. update_distributed_cp complete" << KNRM << std::endl;
        
            leg_duration = 
                ((double)distributed_control_points.size() + order) *
                ts.get_knot_interval();
            
            ts.update_timespan(leg_duration);
        
            // Add the distributed control points and add the past order number of 
            // control points from the previous iteration
            altered_distributed_control_points.clear();
            altered_distributed_control_points = concatenate_distributed_cp();
            std::cout << "[main_server]" << KBLU << 
                " 5. concatenate_distributed_cp complete" << KNRM << std::endl;
        
            // Check whether there is a need to clamp the start when time is 0
            if ((global_search_path[0] - start).norm() < 0.0001 &&
                ts.get_duration_from_start_time() > pow(10,6))
            {
                for (int i = 0; i < order; i++)
                    altered_distributed_control_points.insert(
                        altered_distributed_control_points.begin(), start);
            }

            // Clamp the end
            for (int i = 0; i < order; i++)
                altered_distributed_control_points.insert(
                    altered_distributed_control_points.end(), end);
        }

        /** @brief Outdated **/
        // extract_direct_goal_radius();
        // std::cout << "[main_server]" << KBLU << 
        //     " 4. extract_direct_goal complete" << KNRM << std::endl;
        // std::cout << "[main_server] direct_goal: " << KBLU << 
        //     direct_goal.transpose() << KNRM << std::endl;
        
        /** @brief Outdated **/
        // extract_direct_goal_velocity();
        // std::cout << "[main_server] direct_goal: " << KBLU << 
        //     direct_goal.transpose() << KNRM << std::endl;      
        
        // std::cout << "[main_server] get_duration_from_start_time() " << KBLU << 
        //     ts.get_duration_from_start_time() << KNRM << std::endl;

        // Update the [timespan] in trajectory_server node, extending it relative
        // to our queried control point

        // Get the valid control points that we can use
        // vector<Eigen::Vector3d> acceptable_cp =
        //     ts.get_valid_cp_vector(altered_distributed_control_points);

        if (!previous_search_points_time.empty())
        {
            if (ts.get_running_time() >= previous_search_points_time[0])
            {
                if (previous_search_points_time.size() > 1)
                {
                    previous_search_points_time.erase(previous_search_points_time.begin());
                    previous_search_points.erase(previous_search_points.begin());
                }
            } 
        }
        
        ts.update_control_points(altered_distributed_control_points);
        std::cout << "[main_server]" << KBLU << 
            " 6. update_control_points complete" << KNRM << std::endl;

        vector<Eigen::Vector3d> local_points =
            get_bspline_control_points(duration_secs);
        
        if ((int)local_points.size() > 0 && (int)local_points.size() > order*2 &&
            ts.get_duration_from_start_time() > get_bspline_knot_interval() &&
            ts.get_duration_from_start_time() < pow(10,6))
        {
            // std::cout << "[main_server] local_points ";
            // for (int j = 0; j < local_points.size(); j++)
            // {
            //     std::cout << "[" << local_points[j].transpose() << "] "; 
            // }
            // std::cout << std::endl;

            run_optimization(local_points);
            update_bspline_local_control_points(
                duration_secs, optimized_control_points);
            std::cout << "[main_server]" << KBLU << 
                " 7. optimized and update complete" << KNRM << std::endl;
        }
        else
        {
            std::cout << "[main_server]" << 
                KRED << " 7. Not Running Optimization" << KNRM << std::endl;
        }

        altered_distributed_control_points.clear();
        altered_distributed_control_points = ts.get_bspline_control_points();

    }

    void main_server::run_optimization(
        vector<Eigen::Vector3d> local_control_points)
    {
        /** @brief Previous implementation summary **/
        // vector<Eigen::Vector3d> solver(
        //     vector<Eigen::Vector3d> global_spline, 
        //     vector<Eigen::Vector3d> reference_spline, 
        //     vector<double> time_points, 
        //     pcl::PointCloud<pcl::PointXYZ>::Ptr obs, 
        //     double dt, double max_acc,
        //     double protected_zone,
        //     double min_height, double max_height);

        // global_spline = get_bspline_control_points (finite)
        // reference_spline = fixed (finite)
        // time_points = knots
        // obs = local_cloud
        // dt = knot_interval
        // max_acc = maximum acceleration
        // protected_zone = obs_threshold
        // min_height = min_height
        // max_height = max_height

        // Seems like LBFGS-B is using "More Thuente" for line search by default

        // Set up parameters
        LBFGSBParam<double> param = ou.setup_lbfgs_param();
                
        vector<double> time_points = 
            get_bspline_knots(duration_secs);

        int number_of_col = local_control_points.size();
        int number_of_row = 3; // Should be 3 for XYZ

        // We need to compress g_spline to a single array
        // g_spline -> g_single_array
        // MatrixXd to vector<double>
        // col*0 + (0 to col-1) = x 
        // col*1 + (0 to col-1) = y
        // col*2 + (0 to col-1) = z
        // Initial guess = g_single_array which is x
        Eigen::VectorXd x = 
            VectorXd::Zero(number_of_col * number_of_row);
        x = ou.vect_vector3_to_single(local_control_points);

        // Create solver and function object
        LBFGSBSolver<double> solver(param);  // New solver class
        trajectory_server::optimizer opt(number_of_col * number_of_row);

        // Bounds
        // Setup Lower Bound and upper bound for now we take z to be clamped
        Eigen::VectorXd lb = VectorXd::Zero(number_of_col * number_of_row);
        Eigen::VectorXd ub = VectorXd::Zero(number_of_col * number_of_row);
        for (int i = 0; i < number_of_row; i++)
        {
            for (int j = 0; j < number_of_col; j++)
            {
                // if (i != 2)
                // {
                    // If representing x and y
                    // Load in max and min double lower and upper bound
                    lb(number_of_col*i + j) = -DBL_MAX;
                    ub(number_of_col*i + j) = DBL_MAX;
                // }
                // else
                // {
                //     // If representing z
                //     // Load in Lower and Upper bound
                //     // lb(number_of_col*i + j) = min_height;
                //     // ub(number_of_col*i + j) = max_height;
                //     lb(number_of_col*i + j) = -DBL_MAX;
                //     ub(number_of_col*i + j) = DBL_MAX;
                // }

            }
        }

        // x will be overwritten to be the best point found
        double fx = 0; // Cost

        std::lock_guard<std::mutex> others_lock(other_traj_mutex);
        opt.set_params_and_data(
            _uav_idx,
            get_bspline_knot_interval(), _max_acc, 
            obs_threshold, local_control_points, 
            local_cloud, time_points, other_agents,
            get_bspline_chronos_start_time());
        opt.set_weights(_weight_vector);

        std::cout << "[main_server] single array size " << KBLU << 
            x.size() << KNRM << std::endl;
        
        
        int iter = solver.minimize(opt, x, fx, lb, ub);

        // single_array -- back into --> global_spline
        // vector<double> to vector<Vector3d>

        std::cout << "[main_server] Iterations " << 
            KBLU << iter << KNRM << "! F(x) " << 
            KBLU << fx << KNRM << "!" << std::endl;
        
        optimized_control_points.clear();
        bool early_break = false;
        std::cout << "[main_server] optimized_control_points ";
        for (int j = 0; j < number_of_col; j++)
        {
            // Squash the nan bug that will kill the program if no solution is found
            if (isnan(x[number_of_col * 0 + j]) || 
                isnan(x[number_of_col * 1 + j]) ||
                isnan(x[number_of_col * 2 + j]))
            {
                std::cout << "[main_server] nan caught on " << 
                    KRED << j << KNRM << std::endl;
                early_break = true;
                break;
            }
            // Add the vector3d to the global_spline vector
            Eigen::Vector3d single_point = Eigen::Vector3d(
                x[number_of_col * 0 + j], x[number_of_col * 1 + j], x[number_of_col * 2 + j]);
            
            // std::cout << "[" << single_point.transpose() << "] "; 
            optimized_control_points.push_back(single_point);
        }
        std::cout << std::endl;

        if (early_break)
        {
            std::cout << "[main_server] " << 
                KRED << "Reject Optimization" << KNRM << std::endl;
            optimized_control_points = local_control_points;
        }

    }

    void main_server::update_other_agents(
        int idx, vector<Eigen::Vector3d> cp, vector<double> knots,
        system_clock::time_point origin_time)
    {
        std::lock_guard<std::mutex> others_lock(other_traj_mutex);
        if (!other_agents.empty())
        {
            int vector_index = -1;
            for (int i = 0; i < other_agents.size(); i++)
            {
                if (idx == other_agents[i].id)
                {
                    vector_index = i;
                    break;
                }
            }

            if (vector_index < 0)
            {
                trajectory_server::
                    optimization_utils::other_agents_traj agent;

                agent.id = idx;
                agent.cp = cp;
                agent.knots = knots;
                agent.origin = origin_time;

                other_agents.push_back(agent);
                return;
            }

            
            other_agents[vector_index].cp = cp;
            other_agents[vector_index].knots = knots;
            other_agents[vector_index].origin = origin_time;
            return;
        }
        else
        {
            trajectory_server::
                optimization_utils::other_agents_traj agent;

            agent.id = idx;
            agent.cp = cp;
            agent.knots = knots;
            agent.origin = origin_time;

            other_agents.push_back(agent);
            return;
        }
    }

    /** @brief Outdated **/
    /** @brief Extract the direct goal for the planner that is according to maximum velocity */
    // void main_server::extract_direct_goal_velocity()
    // {
    //     previous_search_points.clear();
        
    //     // Get the distance between the previous_search_points 
    //     // and find the new time vector
    //     double local_distance_permissable =
    //         max_vel * duration_secs;
    //     std::cout << "[main_server] local_distance_permissable: " << KGRN << 
    //         local_distance_permissable << KNRM << std::endl;
        
    //     vector<Eigen::Vector3d> temporary_vector = global_search_path;

    //     double total_distance = 0.0, distance_left = 0.0, distance;
        
    //     // Find the intersection between any of the legs and the sphere perimeter
    //     intersection_idx = -1;
    //     for (int i = 0; i < temporary_vector.size(); i++)
    //     {
    //         // Find the distance between the 2 points
    //         distance = 
    //             (temporary_vector[i+1] - temporary_vector[i]).norm();
    //         std::cout << "[main_server] " << i << " distance: " << KGRN << 
    //             distance << KNRM << std::endl;
    //         // Check whether the next distance is larger than the local_distance
    //         if (total_distance + distance > local_distance_permissable)
    //         {
    //             // Find out the remaining distance
    //             distance_left = local_distance_permissable - total_distance;
    //             // Intersection index same as the convention used in extract_direct_goal_radius()
    //             intersection_idx = i+1;
    //             std::cout << "[main_server] distance_left: " << KGRN << 
    //                 distance_left << KNRM << std::endl;
    //             // Find the direct goal from the distance left and the vector
    //             direct_goal = ((temporary_vector[i+1] - temporary_vector[i]) / 
    //                 (temporary_vector[i+1] - temporary_vector[i]).norm()) * 
    //                 distance_left + temporary_vector[i];
    //             break;
    //         }

    //         // For the next iteration. add the distance to the total distance
    //         total_distance += distance;
    //     }
    //     std::cout << "[main_server] total_distance: " << KGRN << 
    //         total_distance << KNRM << std::endl;

    //     // If the goal is reached, the distance is smaller than the permissable distance
    //     if (total_distance < local_distance_permissable && 
    //         intersection_idx < 0)
    //     {
    //         direct_goal = end;
    //         intersection_idx = temporary_vector.size()-1;
    //     }
        
    //     // Save previous global_search_path into previous_search_points
    //     for (int i = 0; i < (int)global_search_path.size(); i++)
    //     {
    //         // Do not add the start point
    //         if (i == 0)
    //             continue;

    //         previous_search_points.push_back(global_search_path[i]);
    //     }

    // }

    /** @brief Outdated **/
    /** @brief Extract the direct goal for the planner that is within the search sphere */
    // void main_server::extract_direct_goal_radius()
    // {
    //     double radii;
    //     radii = search_radius;
    //     previous_search_points.clear();

    //     // Find the intersection between any of the legs and the sphere perimeter
    //     intersection_idx = -1;
    //     for (int i = 0; i < (int)global_search_path.size()-1; i++)
    //     {
    //         // outside : previous_search_points[i]
    //         // inside : previous_search_points[i-1]
    //         if (inside_sphere_check(global_search_path[i], current_control_point, search_radius) &&
    //             !inside_sphere_check(global_search_path[i+1], current_control_point, search_radius))
    //         {
    //             intersection_idx = i+1;
    //             std::cout << "[main_server]" << KBLU <<
    //                 " intersection found at idx: " << intersection_idx << 
    //                 KNRM << std::endl;
    //             break;
    //         }
    //     }

    //     /** @brief Print intersection points with sphere */
    //     // std::cout << "start and end global_search_path check [" << intersection_idx-1 <<
    //     //     " and " << intersection_idx << "] out of " <<
    //     //     global_search_path.size() << std::endl;

    //     // We have reached the goal if failed counts is same as global_search_path size
    //     if (intersection_idx < 0)
    //     {
    //         direct_goal = end;
    //         for (int i = 1; i < (int)global_search_path.size(); i++)
    //         {
    //             previous_search_points.push_back(global_search_path[i]);
    //         }
    //         previous_search_points.push_back(end);
    //         return;
    //     }

    //     // direction is inside to outside
    //     Eigen::Vector3d vect1 = global_search_path[intersection_idx] - global_search_path[intersection_idx-1];
    //     Eigen::Vector3d vect1_direction = vect1 / vect1.norm();
    //     // direction is center to inside
    //     Eigen::Vector3d vect2 = global_search_path[intersection_idx-1] - current_control_point;
    //     double a = pow(vect1.x(), 2) + pow(vect1.y(), 2) + pow(vect1.z(), 2);
    //     double b = 2 * (vect1.x() * vect2.x() + 
    //         vect1.y() * vect2.y() + 
    //         vect1.z() * vect2.z());
    //     double c = pow(current_control_point.x(), 2) +
    //         pow(current_control_point.y(), 2) +
    //         pow(current_control_point.z(), 2) +
    //         pow(global_search_path[intersection_idx-1].x(), 2) +
    //         pow(global_search_path[intersection_idx-1].y(), 2) +
    //         pow(global_search_path[intersection_idx-1].z(), 2) -
    //         2 * (global_search_path[intersection_idx-1].x() * current_control_point.x() +
    //         global_search_path[intersection_idx-1].y() * current_control_point.y() +
    //         global_search_path[intersection_idx-1].z() * current_control_point.z()) -
    //         pow(radii,2);

    //     double u_p = (-b + sqrt(pow(b,2) - 4 * a * c)) / (2 * a);

    //     Eigen::Vector3d query1 = vect1 * u_p + global_search_path[intersection_idx-1];

    //     direct_goal = query1;

    //     // Save previous global_search_path into previous_search_points
    //     for (int i = 0; i < (int)global_search_path.size(); i++)
    //     {
    //         /** @brief Push back the intersection points as the final points */
    //         // if (intersection_idx-1 == i)
    //         // {
    //         //     previous_search_points.push_back(query1);
    //         //     break;
    //         // }

    //         // Do not add the start point
    //         if (i == 0)
    //             continue;

    //         previous_search_points.push_back(global_search_path[i]);
    //     }

    //     return;
    // }
}