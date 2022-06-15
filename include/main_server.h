/*
* main_server.h
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

#ifndef MAIN_SERVER_H
#define MAIN_SERVER_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <unistd.h>
#include <cstdlib>
#include <thread>   
#include <mutex>
#include <cfloat>
#include <cmath>
#include <Eigen/Dense>
#include <random>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include "trajectory_server.h"
#include "rrt_server_module.h"
#include "optimization_server.h"

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace std::chrono;
using namespace std;
using namespace Eigen;
using namespace trajectory;
using namespace rrt_server;
using namespace LBFGSpp;

namespace trajectory_server
{
    class main_server
    {
        private:

            /** @brief Front end RRT Server */
            rrt_server::rrt_server_node fe_rrt_server;
            rrt_server::rrt_utility ru; 
            trajectory_server::bspline_server ts;
            trajectory_server::optimization_utils ou;

            vector<trajectory_server::
                optimization_utils::other_agents_traj> other_agents;

            Eigen::Vector3d start, end;
            Eigen::Vector3d direct_goal, current_control_point;
            Eigen::Vector3d current_pose;

            vector<double> _weight_vector;

            vector<Eigen::Vector3d> previous_search_points;
            vector<Eigen::Vector3d> global_search_path;

            vector<Eigen::Vector4d> no_fly_zone;

            /** @brief Start and end time span to feed into bspline */
            vector<Eigen::Vector3d> time_span; 
            vector<Eigen::Vector3d> distributed_control_points, altered_distributed_control_points;
            vector<Eigen::Vector3d> control_points, optimized_control_points;
            vector<double> previous_search_points_time;

            int order, duration_secs, _uav_idx;

            double min_height, max_height;
            double obs_threshold;
            double search_radius;
            double intersection_idx;
            double max_vel, _max_acc;
            double sub_runtime_error, runtime_error;
            double leg_duration;

            std::mutex search_points_mutex;
            std::mutex control_points_mutex;
            std::mutex pose_update_mutex;
            std::mutex cloud_mutex;
            std::mutex other_traj_mutex;

            bool bs_timer_started = false;

            std::random_device dev;

            pcl::PointCloud<pcl::PointXYZ>::Ptr local_cloud;

        public:
            main_server()
            {
                local_cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(
                    new pcl::PointCloud<pcl::PointXYZ>());
            }

            ~main_server(){}

            /** @brief Outdated **/
            /** @brief Extract the direct goal for the planner that is according to maximum velocity */
            // void extract_direct_goal_velocity();

            /** @brief Outdated **/
            /** @brief Extract the direct goal for the planner that is within the search sphere */
            // void extract_direct_goal_radius();

            /** @brief Use this function wisely, since check and update may cause an infinite loop
            * If a bad data is given to the rrt node and it cannot complete the validity check
            */
            void check_and_update_search(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs,  Eigen::Vector3d first_cp);

            /** @brief Construct the search path from RRT search and from its shortened path */
            void generate_search_path(pcl::PointCloud<pcl::PointXYZ>::Ptr obs);

            /** @brief Get the local control points from the RRT module (via distribution) */
            vector<Eigen::Vector3d> get_local_control_points();

            /** @brief Run the whole algorithm to acquire the control points */
            void complete_path_generation();

            void run_optimization(
                vector<Eigen::Vector3d> local_control_points);

            void update_other_agents(
                int idx, vector<Eigen::Vector3d> cp, vector<double> knots,
                system_clock::time_point origin_time);

            /** @brief Update the local cloud data */
            void set_local_cloud(
                pcl::PointCloud<pcl::PointXYZ>::Ptr _local_cloud)
            {
                std::lock_guard<std::mutex> cloud_lock(cloud_mutex);

                if (!local_cloud->empty())
                    local_cloud->points.clear();
                local_cloud = _local_cloud;
            }

            /** @brief Initialize the Bspline server */
            void initialize_bspline_server(
                int _order, double _duration_secs, 
                double _command_interval, int _knot_div,
                double _max_vel)
            {
                order = _order;
                duration_secs = _duration_secs;
                ts.init_bspline_server(_order, _duration_secs, _command_interval, _knot_div);
                max_vel = _max_vel;
            }

            /** @brief Initialize the RRT server */
            void initialize_rrt_server(
                double _sub_runtime_error, double _runtime_error,
                double xybuffer, double zbuffer, double passage_size)
            {
                fe_rrt_server.setup_buffers(xybuffer, zbuffer, passage_size);

                sub_runtime_error = _sub_runtime_error;
                runtime_error = _runtime_error;
            }

            /** @brief Initialize the optimization server */
            void initialize_opt_server(
                vector<double> weight_vector, double max_acc, int uav_idx)
            {
                _uav_idx = uav_idx;
                if (weight_vector.size() != 5)
                    return;
                // Index for weights
                // 0. _weight_smooth = weight_vector[0]; 
                // 1. _weight_feas = weight_vector[1];
                // 2. _weight_term = weight_vector[2];
                // 3. _weight_static = weight_vector[3];
                // 4. _weight_reci = weight_vector[4];
                _weight_vector = weight_vector;
                _max_acc = max_acc;
            }

            /** @brief Restart and reinitialize the start and end goals */
            void reset_goal_points(
                Eigen::Vector3d s, Eigen::Vector3d e, 
                double protected_zone, double search_radii)
            {
                previous_search_points.clear();
                time_span.clear();
                control_points.clear();
                global_search_path.clear();
                distributed_control_points.clear();
                start = s; end = e; 
                obs_threshold = protected_zone;
                search_radius = search_radii;
            }

            /** @brief Concatenate the control points (overlapping + distributed points) */
            vector<Eigen::Vector3d> concatenate_distributed_cp()
            {
                vector<Eigen::Vector3d> new_cp_vector;

                new_cp_vector = ts.get_overlapping_cp();
                for (int i = 0; i < distributed_control_points.size(); i++)
                    new_cp_vector.push_back(distributed_control_points[i]);
                
                return new_cp_vector;
            }

            /** @brief Check whether the point is inside the sphere */
            bool inside_sphere_check(Eigen::Vector3d point, 
                Eigen::Vector3d sphere_center, double radii)
            {
                // get the distance vector
                Eigen::Vector3d dv = sphere_center - point;
                double a = 1.0; 
                double b = 1.0; 
                double inv_a2 = 1 / a / a; 
                double inv_b2 = 1 / b / b; 
                // get the ellipsoidal distance 
                double e_d = sqrt(pow(dv.z(),2) * inv_a2 + 
                    (pow(dv.x(),2) + pow(dv.y(),2)) * inv_b2); 
                // get the distance error
                double d_e = radii - e_d; 
                if (d_e > 0)
                    return true;
                else 
                    return false;
            }

            /** @brief Return the distributed control points */
            vector<Eigen::Vector3d> get_distributed_control_points()
            {
                return distributed_control_points;
            }

            /** @brief Get the command from the Bspline */
            bspline_server::pva_cmd update_bs_path_get_command() 
            {
                return ts.update_get_command_by_time();
            }

            /** @brief Get the end time from the bspline server **/
            double get_end_time()
            {
                return ts.get_end_time_for_path();
            }

            /** @brief Get the start time from the bspline server **/
            double get_start_time()
            {
                return ts.get_start_time_for_path();
            }

            /** @brief Reset the start time for the bspline server **/
            void reset_bspline_start_time()
            {
                ts.reset_start_time();
            }

            /** @brief Start the Bspline timer */
            void start_module_timer()
            {
                ts.start_bspline_time();
            }

            /** @brief Get time from bspline module */
            double get_bspline_time()
            {
                return ts.get_running_time();
            }

            double get_bspline_knot_interval()
            {
                return ts.get_knot_interval();
            }

            time_point<std::chrono::system_clock> 
                get_bspline_chronos_start_time()
            {
                return ts.get_chronos_start_time();
            }

            vector<Eigen::Vector3d> get_bspline_control_points(double time)
            {
                return ts.get_finite_bs_control_points(time);
            }

            vector<double> get_bspline_knots(double time)
            {
                return ts.get_current_knots(time);
            }

            void update_bspline_local_control_points(
                double time, vector<Eigen::Vector3d> control_points)
            {
                ts.update_local_control_points(
                    time, control_points);
            }
    };

};

#endif