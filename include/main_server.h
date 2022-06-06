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

namespace trajectory_server
{
    class main_server
    {
        private:
            /** @brief Front end RRT Server */
            rrt_server::rrt_server_node fe_rrt_server;
            trajectory_server::bspline_server ts;
            rrt_server::rrt_utility ru; 

            Eigen::Vector3d start, end;
            Eigen::Vector3d direct_goal, current_control_point;
            Eigen::Vector3d current_pose;

            vector<Eigen::Vector3d> previous_search_points;
            vector<Eigen::Vector3d> global_search_path;
            vector<Eigen::Vector3d> local_search_path;

            vector<Eigen::Vector4d> no_fly_zone;

            /** @brief Start and end time span to feed into bspline */
            vector<Eigen::Vector3d> time_span; 
            vector<Eigen::Vector3d> distributed_control_points;
            vector<Eigen::Vector3d> control_points, optimized_control_points;

            double min_height, max_height;
            double obs_threshold;
            double search_radius;
            double intersection_idx;
            double max_vel;
            double sub_runtime_error, runtime_error;

            std::mutex search_points_mutex;
            std::mutex control_points_mutex;
            std::mutex pose_update_mutex;
            std::mutex cloud_mutex;

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

            /** @brief Extract the direct goal for the planner that is within the search sphere */
            void extract_direct_goal();

            /** @brief Use this function wisely, since check and update may cause an infinite loop
            * If a bad data is given to the rrt node and it cannot complete the validity check
            */
            void check_and_update_search(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs,  Eigen::Vector3d first_cp);

            /** @brief Construct the search path from RRT search and from its shortened path */
            void generate_search_path(pcl::PointCloud<pcl::PointXYZ>::Ptr obs);

            void complete_path_generation();

            void set_local_cloud(
                pcl::PointCloud<pcl::PointXYZ>::Ptr _local_cloud)
            {
                std::lock_guard<std::mutex> cloud_lock(cloud_mutex);

                if (!local_cloud->empty())
                    local_cloud->points.clear();
                local_cloud = _local_cloud;
            }

            void initialize_bspline_server(
                int _order, double _duration_secs, 
                double _command_interval, int _knot_div,
                double _max_vel)
            {
                ts.init_bspline_server(_order, _duration_secs, _command_interval, _knot_div);
                max_vel = _max_vel;
            }

            void initialize_rrt_server(
                double _sub_runtime_error, double _runtime_error)
            {
                sub_runtime_error = _sub_runtime_error;
                runtime_error = _runtime_error;
            }

            void reset_goal_points(
                Eigen::Vector3d s, Eigen::Vector3d e, 
                double protected_zone, double search_radii)
            {
                previous_search_points.clear();
                time_span.clear();
                control_points.clear();
                start = s; end = e; 
                obs_threshold = protected_zone;
                search_radius = search_radii;
            }

            vector<Eigen::Vector3d> get_local_control_points(double max_vel) 
            {
                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);
                local_search_path.clear();
                for (int i = 0; i < (int)global_search_path.size(); i++)
                {
                    /** @brief Push back the intersection points as the final points */
                    if (intersection_idx-1 == i)
                    {
                        local_search_path.push_back(direct_goal);
                        break;
                    }
                    local_search_path.push_back(previous_search_points[i]);
                }

                return ts.get_redistributed_cp_vector(
                    current_control_point,
                    local_search_path, max_vel);
            }

            bspline_server::pva_cmd update_distributed_cp(double max_vel) 
            {
                distributed_control_points.clear();
                distributed_control_points = get_local_control_points(max_vel);
            }

            bspline_server::pva_cmd update_bs_path_get_command(double max_vel) 
            {
                trajectory_server::bspline_server::pva_cmd cmd_by_time;

                vector<Eigen::Vector3d> acceptable_cp =
                    ts.get_valid_cp_vector(distributed_control_points);

                return cmd_by_time;
            }

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

            vector<Eigen::Vector3d> get_distributed_control_points()
            {
                return distributed_control_points;
            }
    };

};

#endif