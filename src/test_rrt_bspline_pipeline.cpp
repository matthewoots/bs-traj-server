/*
* test_rrt_bspline_pipeline.cpp
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
#include "main_server.h"
#include <random>
// #include <iomanip>

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
    trajectory_server::main_server ms;

    std::random_device dev;
    std::mt19937 generator(dev());
    std::uniform_real_distribution<double> dis_middle(-1.0, 1.0);
    std::uniform_real_distribution<double> dis_normal(0.0, 1.0);

    Eigen::Vector3d start = Eigen::Vector3d(-5,5,1.5);
    Eigen::Vector3d goal = Eigen::Vector3d(5,-5,1.9);

    double obs_threshold = 0.2;
    double search_radius = 9.0;

    double _sub_runtime_error = 0.02;
    double _runtime_error = 0.12;

    double _min_height = 1.3;
    double _max_height = 2.0;

    double _xybuffer = 1.0;
    double _zbuffer = 0.5;
    double _passage_size = 12.0;

    double _max_vel = 3.0;

    double _duration = 8.0;
    double _cmd_update_hz = 10;
    double _cmd_update_interval = 1/_cmd_update_hz;

    double path_generation_interval = 0.8;
    double time_to_add_cloud_iter = 2;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(
        new pcl::PointCloud<pcl::PointXYZ>());
    // size of cloud = 1 x 1 x 1 block
    double cube_size = 1.0; 
    Eigen::Vector3d middle = (start + goal) / 2;
    double offset = - cube_size / 2;
    double resolution = obs_threshold / 2;
    int index = (int)(1 / resolution);


    for (int i = 0; i < index; i++)
    {
        for (int j = 0; j < index; j++)
        {
            for (int k = 0; k < index; k++)
            {
                pcl::PointXYZ point;
                point.x = middle.x() + offset + i*resolution;
                point.y = middle.y() + offset + j*resolution;
                point.z = middle.z() + offset + k*resolution;
                cloud->points.push_back(point);
            }
        }
    } 

    time_point<std::chrono::system_clock> test_cycle_start = system_clock::now();   

    std::cout << "[pipeline] " << KBLU << 
        "reset_goal_points" << KNRM << std::endl;
    ms.reset_goal_points(
        start, goal, obs_threshold, search_radius);
    // ms.initialize_bspline_server(
    //  _order, _traj_duration_secs, 1/_cmd_update_hz, 
    //  _des_knot_div, _max_vel);
    std::cout << "[pipeline] " << KBLU << 
        "initialize_bspline_server" << KNRM << std::endl;
    ms.initialize_bspline_server(
        5, _duration, _cmd_update_interval, 5, _max_vel);
    std::cout << "[pipeline] " << KBLU << 
        "initialize_rrt_server" << KNRM << std::endl;
    ms.initialize_rrt_server(
        _sub_runtime_error, _runtime_error,
        _xybuffer, _zbuffer, _passage_size,
        _min_height, _max_height);

    std::cout << "[pipeline] " << KBLU << 
        "complete_path_generation" << KNRM << std::endl;
    ms.complete_path_generation();

    std::cout << "[pipeline] " << KBLU << 
        "start_module_timer" << KNRM << std::endl;
    ms.start_module_timer();

    int time_step_ms = (int)floor(_cmd_update_interval * 1000);
    double time = 0.0;
    std::cout << "[pipeline] duration: " << KBLU << 
        ms.get_end_time() << KNRM << std::endl;
    
    // i represents the command timer, j represents the path_generation timer
    int i = 1, j = 1;
    std::this_thread::sleep_for(std::chrono::milliseconds(time_step_ms));
    while (time < ms.get_end_time() - _cmd_update_interval)
    {
        time = i * _cmd_update_interval;


        if (time - j * path_generation_interval >= 0)
        {
            j++;
            if (j == time_to_add_cloud_iter)
            {
                std::cout << "[pipeline] Added cloud" << std::endl;
                ms.set_local_cloud(cloud);
            }
            
            std::cout << std::endl;
            
            std::cout << KRED << 
                "[pipeline] Starting complete_path_generation" << 
                KNRM << std::endl;
            time_point<std::chrono::system_clock> test_cycle_start = system_clock::now();

            ms.complete_path_generation();
            
            auto test_time_diff = duration<double>(system_clock::now() - test_cycle_start).count();
            std::cout << KRED << 
                "[pipeline] Completed complete_path_generation" << 
                KNRM << " in " << KBLU << test_time_diff*1000 << KNRM <<
                "ms" << std::endl;

            std::cout << std::endl;          
        }

        trajectory_server::bspline_server::pva_cmd cmd;
        cmd = ms.update_bs_path_get_command();

        // IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
        // std::cout << "[pipeline] time: " << 
        //     KRED << std::setprecision(3) << ms.get_start_time()  << KNRM <<
        //     " " << 
        //     KBLU <<  std::setprecision(3) << cmd.t  << KNRM <<
        //     " " << 
        //     KRED << std::setprecision(3) << ms.get_end_time()  << KNRM <<
        //     " position: " << 
        //     KBLU << cmd.p.transpose().format(CleanFmt) << KNRM << 
        //     " total velocity: " << 
        //     KBLU << std::setprecision(3) << cmd.v.norm() << KNRM << std::endl;
        
        double offset_ms = 0.0;
        if (ms.get_bspline_time() - time > 0)
            offset_ms = (ms.get_bspline_time() - time) * 1000;
        std::this_thread::sleep_for(std::chrono::milliseconds(time_step_ms - (int)offset_ms));
        i++;
    }

    // we have to close the trajectory everytime we finish
    
    return 1;
}