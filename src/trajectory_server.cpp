/*
* trajectory_server.cpp
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

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

namespace trajectory_server
{
    bspline_trajectory::bs_pva_state_3d bspline_server::get_bs_path
        (vector<Eigen::Vector3d> cp)
    {
        bspline_trajectory::bs_pva_state_3d pva3;
        int acceptable_cp_size = bsu.get_valid_cp_size(
            knot_interval, order, timespan);
        
        std::cout << "[tserver]" <<
            " cp_size/acceptable " <<
            cp.size() << "/" << acceptable_cp_size << std::endl;
        
        // Crop the control points to what can be given to the bspline
        vector<Eigen::Vector3d> acceptable_cp;
        for (int i = 0; i < acceptable_cp_size; i++)
            acceptable_cp.push_back(cp[i]);

        pva3 = bsu.get_uni_bspline_3d(
            order, timespan, acceptable_cp, knot_div);

        std::cout << "[tserver]" <<
            " path_size: " <<
            pva3.pos.size() << std::endl;

        return pva3;
    }
    
    bspline_server::pva_cmd bspline_server::get_command_on_path()
    {
        bspline_server::pva_cmd pva;

        return pva;
    }

    bool bspline_server::valid_cp_count_check(size_t cp_size)
    {
        return ((int)cp_size > (knot_size + order));
    }
    
}