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
using namespace std;

namespace trajectory_server
{
    class bspline_server
    {
        public:
            bspline_server(int _order, double duration) 
            {
                order = _order;
                stime = std::chrono::system_clock::now();
            }
            ~bspline_server() {std::cout << "close bs_server" << std::endl;}

        private:
            /** @brief Parameters for libbspline */ 
            int order;
            vector<double> timespan;
            std::chrono::time_point<std::chrono::system_clock> stime; // start time for bspline server in time_t

    };

};

#endif