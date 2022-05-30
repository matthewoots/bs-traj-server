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
#include <Eigen/Dense>
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

            Eigen::Vector3d start, end;
            Eigen::Vector3d direct_goal;
            Eigen::Vector3d current_pose;

            /** @brief Start and end time span to feed into bspline */
            vector<Eigen::Vector3d> time_span; 
            vector<Eigen::Vector3d> previous_search_points;
            vector<Eigen::Vector4d> no_fly_zone;
            vector<Eigen::Vector3d> control_points;

            double min_height, max_height;

            std::mutex search_points_mutex;
            std::mutex control_points_mutex;
            std::mutex pose_update_mutex;

        public:
            main_server(){}

            ~main_server(){}

            void reset_goal_points(
                Eigen::Vector3d s, Eigen::Vector3d e)
            {
                previous_search_points.clear();
                time_span.clear();
                control_points.clear();
                start = s; end = e; 
            }

            void check_and_update_search()
            {
                std::lock_guard<std::mutex> pose_lock(pose_update_mutex);
                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);

                double nearest_distance = DBL_MAX;
                int idx = -1; 
                for (int i = 0; i < previous_search_points.size(); i++)
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
                check_and_update_search();

                vector<Eigen::Vector3d> path = fe_rrt_server.find_rrt_path(
                    previous_search_points, obs_pcl, start, end, 
                    no_fly_zone, min_height, max_height);

                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);
                if (path.empty())
                    return;
                
                previous_search_points.clear();
                previous_search_points = path;

                return;
            }

            bspline_server::pva_cmd update_bs_path_get_command()
            {
                std::lock_guard<std::mutex> search_points_lock(search_points_mutex);

                trajectory_server::bspline_server::pva_cmd cmd_by_time;

                return cmd_by_time;
            }

    };

};

#endif