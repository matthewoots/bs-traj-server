/*
* trajectory_server.h
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

#ifndef TRAJECTORY_SERVER_H
#define TRAJECTORY_SERVER_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <unistd.h>
#include <cstdlib>
#include <thread>   
#include <mutex>
#include <Eigen/Dense>

#include "bspline_utils.hpp"

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace Eigen;
using namespace trajectory;
using namespace std::chrono;
using namespace std;

namespace trajectory_server
{
    class test_utils
    {
        public:
        void test_chronos()
        {
            time_point<std::chrono::system_clock> start = system_clock::now();
            std::this_thread::sleep_for(std::chrono::milliseconds(1)); // sleep for milliseconds
            time_point<std::chrono::system_clock> end = system_clock::now();
            time_t end_time_t = system_clock::to_time_t(end);
            auto test_time_diff = duration<double>(end - start).count();
            std::cout << "[tserver]" << " test chronos: " << 
                KGRN << test_time_diff << KNRM << "s" << std::endl;
        }

        void print_time(time_point<std::chrono::system_clock> query, string topic)
        {
            time_point<std::chrono::system_clock> end = system_clock::now();
            auto test_time_diff = duration<double>(end - query).count();
            std::cout << topic << " " << 
                KGRN << test_time_diff << KNRM << "s" << std::endl;
        }
    };

    class bspline_server
    {
        public:

            struct pva_cmd
            {
                double t; // Relative to start time 
                Eigen::Vector3d p;
                Eigen::Vector3d v;
                Eigen::Vector3d a;
                double yaw;
            };

            trajectory_server::test_utils tu;
            bspline_trajectory bsu;
            common_trajectory_tool ctt;
        
            bspline_server(){}
            
            ~bspline_server() 
            {
                std::cout << "[tserver]" << 
                    KRED << " close bs_server" << KNRM << std::endl;
            }

            /** @brief Initialize the parameters for the Bspline server
             * without starting the timer
            **/
            void init_bspline_server(
                int _order, double _duration_secs, double _command_interval, int _knot_div);


            /** @brief Start the Bspline timer and saves the start time **/
            void start_bspline_time();


            /** @brief outdated **/
            /** @brief Get the corrected control points that fits the knot span
             * @param cp raw control point data that has not yet been cropped to the
             * acceptable size
            **/
            // vector<Eigen::Vector3d> get_valid_cp_vector(vector<Eigen::Vector3d> cp);


            /** @brief outdated **/
            /** @brief get_bs_path, update_get_command_on_path_by_idx and update_bs_path
             * These 3 functions comes together, this is because update_bs_path have to generate
             * the bspline path and the command functions has to search for the idx
            */
            // bspline_trajectory::bs_pva_state_3d get_bs_path(
            //     vector<Eigen::Vector3d> cp);
            // void update_bs_path(vector<Eigen::Vector3d> cp);
            // bspline_server::pva_cmd update_get_command_on_path_by_idx();


            /** @brief update_get_command_on_path_by_time
             * Do not have to generate the whole path again, just use the matrix to find the command
             * at that point in time
            */
            bspline_server::pva_cmd update_get_command_by_time();


            /** @brief Get the uniform distributed control point vector from 
             * RRT points and the maximum velocity **/
            vector<Eigen::Vector3d> get_redistributed_cp_vector(
                Eigen::Vector3d current_target_cp,
                vector<Eigen::Vector3d> cp, double max_vel);

            /** @brief Return the query control point and update the overlap control 
             * point vectors **/
            Eigen::Vector3d get_current_cp_and_overlap();


            bool valid_cp_count_check(size_t cp_size)
            {
                return ((int)cp_size > (knot_size + order));
            }

            double get_running_time()
            {
                time_point<std::chrono::system_clock> now = system_clock::now();
                return duration<double>(now - stime).count();
            }

            double get_duration_from_start_time()
            {
                return duration<double>(system_clock::now() - stime).count();
            }

            void reset_start_time()
            {
                stime = default_time;
            }

            double get_end_time_for_path()
            {
                return timespan[1];
            }

            double get_start_time_for_path()
            {
                return timespan[0];
            }

            double get_knot_interval()
            {
                return knot_interval;
            }
            
            vector<Eigen::Vector3d> get_overlapping_cp()
            {
                return overlapping_control_points;
            }

            void update_local_control_points(
                double time, vector<Eigen::Vector3d> control_points)
            {

                int valid_count = (int)floor(time / knot_interval);
                int counter = 0;
                int base_idx = current_cp_idx - 1;
                for (int i = current_cp_idx - 1; i < bs_control_points.size(); i++)
                {
                    counter++;
                    if (counter > valid_count)
                        break;
                    bs_control_points[i] = control_points[i - base_idx];
                }
            }

            vector<double> get_current_knots(double time)
            {
                vector<double> finite_knots;

                int valid_count = (int)floor(time / knot_interval);
                int counter = 0;
                for (int i = current_cp_idx - 1; i < bs_control_points.size(); i++)
                {
                    counter++;
                    if (counter > valid_count)
                        break;
                    int idx = i - (current_cp_idx - 1);
                    finite_knots.push_back(current_knot_time + idx * knot_interval);
                }
                return finite_knots;
            }

            void update_control_points(vector<Eigen::Vector3d> control_points)
            {
                std::lock_guard<std::mutex> path_lock(bs_path_mutex);
                bs_control_points.clear();
                bs_control_points = control_points;
            }

            void update_timespan(double total_time)
            {
                duration_secs = total_time;
                std::cout << "[bspline_server] duration_secs: " 
                    << KBLU << duration_secs << KNRM << "s" << std::endl;
                if (!timespan.empty())
                {
                    timespan.clear();
                    timespan.push_back(current_knot_time);
                    timespan.push_back(current_knot_time + duration_secs);
                }
            }

            vector<Eigen::Vector3d> get_finite_bs_control_points(double time)
            {
                vector<Eigen::Vector3d> finite_cp;

                int valid_count = (int)floor(time / knot_interval);
                int counter = 0;
                for (int i = current_cp_idx - 1; i < bs_control_points.size(); i++)
                {
                    counter++;
                    if (counter > valid_count)
                        break;
                    finite_cp.push_back(bs_control_points[i]);
                }

                return finite_cp;
            }

        private:

            std::mutex bs_path_mutex;
            std::mutex time_mutex;

            /** @brief Parameters for libbspline */ 
            int order, knot_div, knot_size;
            int acceptable_cp_size;
            double command_interval, duration_secs, knot_interval;
            double current_knot_time;
            double previous_yaw;
            vector<double> timespan;
            vector<double> original_timespan;

            int overlap, current_cp_idx;

            vector<Eigen::Vector3d> bs_control_points;
            vector<Eigen::Vector3d> overlapping_control_points;

            bspline_trajectory::bs_pva_state_3d pva_state;
            bspline_server::pva_cmd pva_cmd_msg;

            time_point<std::chrono::system_clock> default_time;
            time_point<std::chrono::system_clock> stime; // start time for bspline server in time_t
            vector<time_point<std::chrono::system_clock>> time_span_chronos;

    };

};

#endif