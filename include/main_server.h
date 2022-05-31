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

            /** @brief Start and end time span to feed into bspline */
            vector<Eigen::Vector3d> time_span; 
            vector<Eigen::Vector3d> previous_search_points;
            vector<Eigen::Vector4d> no_fly_zone;
            vector<Eigen::Vector3d> control_points;

            double min_height, max_height;
            double obs_threshold;
            double search_radius;

            std::mutex search_points_mutex;
            std::mutex control_points_mutex;
            std::mutex pose_update_mutex;

            bool bs_timer_started = false;

            std::random_device dev;

        public:
            main_server(){}

            ~main_server(){}

            void initialize_bspline_server(
                int _order, double _duration_secs, double _command_interval, int _knot_div)
            {
                ts.init_bspline_server(_order, _duration_secs, _command_interval, _knot_div);
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

            void extract_direct_goal(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs)
            {
                double radii;
                // If end goal is inside the search sphere
                if ((end - current_control_point).norm() < search_radius)
                    radii = (end - current_control_point).norm();
                // If end goal is outside the search sphere
                else
                    radii = search_radius;

                // Find a RRT path that is quick and stretches to the end point
                vector<Eigen::Vector3d> path = fe_rrt_server.find_rrt_path(
                    previous_search_points, obs, current_control_point, end, 
                    no_fly_zone, min_height, max_height, radii, obs_threshold);
                // If path gives an invalid value, execute some form of emergency
                if (path.empty())
                {
                    return;
                }

                // Find the intersection between any of the legs and the sphere perimeter
                int failed_counts = 0;
                int intersection_idx;
                for (int i = 1; i < (int)path.size(); i++)
                {
                    // outside : previous_search_points[i]
                    // inside : previous_search_points[i-1]
                    if (inside_sphere_check(path[i-1], current_control_point, search_radius) &&
                        !inside_sphere_check(path[i], current_control_point, search_radius))
                    {
                        intersection_idx = i;
                        break;
                    }
                    failed_counts++;
                }

                // We have reached the goal if failed counts is same as path size
                if (failed_counts == path.size()-1)
                {
                    direct_goal = end;
                    return;
                }

                // Calculate the direct goal point from the intersection pair
                // Find the intersection with the sphere
                Eigen::Vector3d vect1 = path[intersection_idx] - path[intersection_idx-1];
                double vect1_norm = vect1.norm();
                Eigen::Vector3d vect2 = current_control_point - path[intersection_idx-1];
                double vect2_norm = vect2.norm();
                double dot = vect1.x()*vect2.x() + vect1.y()*vect2.y() + vect1.z()*vect2.z();
                double extension = dot / vect1_norm;
                double direct_distance_from_start;
                if (extension <= 0.00001 && extension >= -0.00001)
                {
                    double nearest_dist = sqrt(pow(vect2_norm,2) - pow(extension,2));
                    double sub_distance = sqrt(pow(radii,2) - pow(nearest_dist,2)); 
                    direct_distance_from_start = sub_distance + extension;
                }
                else
                {
                    direct_distance_from_start = sqrt(pow(radii,2) - pow(extension,2)); 
                }

                direct_goal = direct_distance_from_start * (vect1 / vect1_norm);
                return;
            }

            void check_and_update_search(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs)
            {
                std::lock_guard<std::mutex> pose_lock(pose_update_mutex);
                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);

                if (previous_search_points.empty())
                    return;

                int last_safe_idx = -1;
                int initial_size = (int)previous_search_points.size();
                for (int i = 1; i < initial_size; i++)
                {
                    if (!ru.check_line_validity_with_pcl(
                        previous_search_points[i], previous_search_points[i-1], obs_threshold, obs))
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

                int new_size = (int)previous_search_points.size();
                double nearest_distance = DBL_MAX;
                int idx = -1; 
                for (int i = 0; i < new_size; i++)
                {
                    double distance = (previous_search_points[i]-current_pose).norm();
                    if (distance < nearest_distance)
                    {
                        nearest_distance = distance;
                        idx = i;
                    }
                }

                for (int i = idx-1; i >= 0; i--)
                    previous_search_points.erase(previous_search_points.begin());
            }


            void generate_search_path(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl)
            {
                // check previous search points and create an updated one
                // We moved forward, so find the closes point from current pose
                check_and_update_search(obs_pcl);

                vector<Eigen::Vector3d> path = fe_rrt_server.find_rrt_path(
                    previous_search_points, obs_pcl, current_control_point, direct_goal, 
                    no_fly_zone, min_height, max_height, search_radius / 2, obs_threshold);

                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);
                if (path.empty())
                    return;
                
                previous_search_points.clear();
                previous_search_points = path;

                return;
            }

            vector<Eigen::Vector3d> get_control_points_bs(double max_vel) 
            {
                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);

                return ts.get_redistributed_cp_vector(previous_search_points, max_vel);
            }

            bspline_server::pva_cmd update_bs_path_get_command(double max_vel) 
            {
                trajectory_server::bspline_server::pva_cmd cmd_by_time;

                vector<Eigen::Vector3d> distributed_cp =
                    get_control_points_bs(max_vel);

                vector<Eigen::Vector3d> acceptable_cp =
                    ts.get_valid_cp_vector(distributed_cp);

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
    };

};

#endif