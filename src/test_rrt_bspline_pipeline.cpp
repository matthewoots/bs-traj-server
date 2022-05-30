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
    double map_size = 6.0;
    double _min_height = 0.0;
    double _max_height = 5.0;
    
    int test_cycles = 2;
    double random_multiplier = 5.0;

    for (int i = 0; i < test_cycles; i++)
    {
        time_point<std::chrono::system_clock> test_cycle_start = system_clock::now();

        Eigen::Vector3d start = 
            Eigen::Vector3d(dis_normal(generator) * map_size, dis_normal(generator) * map_size, 
            dis_normal(generator) * map_size);
        Eigen::Vector3d end = 
            Eigen::Vector3d(dis_middle(generator) * map_size, dis_middle(generator) * map_size, 
            dis_normal(generator) * map_size);
    }     

    // we have to close the trajectory everytime we finish
    
    return 1;
}