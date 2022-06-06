/*
* test_random_points.cpp
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
#include <random>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

int main()
{
    trajectory_server::test_utils tu;

    std::random_device dev;
    std::mt19937 generator(dev());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // start our trajectory server
    // bspline_server(int _order, double _duration_secs, double _command_interval) 
    // 25Hz publish rate for command
    trajectory_server::bspline_server ts;
    
    int test_cycles = 2;
    double random_multiplier = 5.0;

    for (int i = 0; i < test_cycles; i++)
    {
        time_point<std::chrono::system_clock> test_cycle_start = system_clock::now();

        ts.init_bspline_server(5, 8.0, 0.04, 10);
        int total_cp_size = (int)round(dis(generator) * 100);

        vector<Eigen::Vector3d> random_cp;
        for (int j = 0; j < total_cp_size; j++)
            random_cp.push_back(Eigen::Vector3d(
                (dis(generator) - 0.5) * random_multiplier, 
                (dis(generator) - 0.5) * random_multiplier, 
                dis(generator) * random_multiplier));
        
        if (!ts.valid_cp_count_check(random_cp.size()))
        {
            std::cout << KRED << "fail valid cp_size check" << KNRM << std::endl;
            continue;
        }

        bspline_trajectory::bs_pva_state_3d state;
        state = ts.get_bs_path(random_cp);
        
        
        // Print out the time taken for the test cycle
        tu.print_time(test_cycle_start, 
            "test (" + std::to_string(i) + ") runtime:");


        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        double displaced_time = ts.get_duration_from_start_time();
        // Print out the time taken for the displacement
        std::cout << "displaced_time: " << 
            KBLU << displaced_time << KNRM << std::endl;

        trajectory_server::bspline_server::pva_cmd cmd_by_idx;
        trajectory_server::bspline_server::pva_cmd cmd_by_time;
        ts.update_bs_path(random_cp);
        vector<Eigen::Vector3d> acceptable_cp =
            ts.get_valid_cp_vector(random_cp);
        
        double test_command_time = ts.get_duration_from_start_time();
        cmd_by_time = ts.update_get_command_by_time();
        
        double test_command_idx = ts.get_duration_from_start_time();
        cmd_by_idx = ts.update_get_command_on_path_by_idx();

        // Print out the output for update_get_command_on_path_by_idx
        std::cout << "update_get_command_on_path_by_idx: \n" << 
            KBLU << cmd_by_idx.p << KNRM << std::endl <<
            "cmd_time: " << cmd_by_idx.t << std::endl <<
            "cmd_difference: " << cmd_by_idx.t - test_command_idx << std::endl;
        // Print out the output for update_get_command_by_time
        std::cout << "update_get_command_by_time: \n" << 
            KBLU << cmd_by_time.p << KNRM << std::endl <<
            "cmd_time: " << cmd_by_time.t << std::endl <<
            "cmd_difference: " << cmd_by_time.t - test_command_time << std::endl;

        std::cout << "\n";
    }     

    // we have to close the trajectory everytime we finish
    ts.~bspline_server();
    
    return 1;
}