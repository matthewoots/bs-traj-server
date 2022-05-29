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
        
            bspline_server(){}

            void init_bspline_server(
                int _order, double _duration_secs, double _command_interval, int _knot_div) 
            {
                tu.test_chronos();
                command_interval = _command_interval;
                knot_div = _knot_div;
                order = _order;
                stime = std::chrono::system_clock::now();

                // Re-adjust duration so that our knots and divisions are matching
                duration_secs = bsu.get_corrected_duration(
                    command_interval, _duration_secs);

                knot_size = bsu.get_knots_size(
                    command_interval, duration_secs, knot_div);

                knot_interval = duration_secs / knot_size;

                // convert duration from secs to milliseconds
                long int duration_msecs = (long int)(_duration_secs * 1000);
                // initialize the very first pair of time span (start and end)
                time_span_chronos.clear();
                time_span_chronos.push_back(stime);
                time_span_chronos.push_back(
                    stime + milliseconds(duration_msecs)); // from the alias std::chrono::seconds
                
                timespan.clear();
                timespan.push_back(0.0);
                timespan.push_back(duration<double>(time_span_chronos[1] - time_span_chronos[0]).count());
            }
            
            ~bspline_server() 
            {
                std::cout << "[tserver]" << 
                    KRED << " close bs_server" << KNRM << std::endl;
            }

            bspline_trajectory::bs_pva_state_3d get_bs_path(
                vector<Eigen::Vector3d> cp);
            
            /** @brief update_get_command_on_path_by_idx and update_bs_path
             * These 2 functions comes together, this is because update_bs_path have to generate
             * the bspline path and the command functions has to search for the idx
            */
            void update_bs_path(vector<Eigen::Vector3d> cp);
            bspline_server::pva_cmd update_get_command_on_path_by_idx();
            
            /** @brief update_get_command_on_path_by_time
             * Do not have to generate the whole path again, just use the matrix to find the command
             * at that point in time
            */
            bspline_server::pva_cmd update_get_command_by_time(
                vector<Eigen::Vector3d> cp);

            bool valid_cp_count_check(size_t cp_size);

            double get_running_time();

            double get_duration_from_start_time(){
                return duration<double>(system_clock::now() - stime).count();}

            vector<Eigen::Vector3d> get_valid_cp_vector(vector<Eigen::Vector3d> cp);

        private:

            std::mutex bs_path_mutex;
            std::mutex cmd_mutex;

            /** @brief Parameters for libbspline */ 
            int order, knot_div, knot_size;
            double command_interval, duration_secs, knot_interval;
            vector<double> timespan;

            bspline_trajectory::bs_pva_state_3d pva_state;
            bspline_server::pva_cmd pva_cmd_msg;

            time_point<std::chrono::system_clock> stime; // start time for bspline server in time_t
            vector<time_point<std::chrono::system_clock>> time_span_chronos;

    };

};

#endif