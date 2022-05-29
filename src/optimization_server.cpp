/*
* optimization_server.cpp
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

#include "optimization_server.h"

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
    // Solver will solve with the given partial global_spline 
    vector<Eigen::Vector3d> optimization_server::solver(
        vector<Eigen::Vector3d> global_spline, 
        vector<Eigen::Vector3d> reference_spline, 
        vector<double> time_points, 
        pcl::PointCloud<pcl::PointXYZ>::Ptr obs, 
        double dt, double max_acc,
        double protected_zone,
        double min_height, double max_height)
    {
        // Seems like LBFGS-B is using "More Thuente" for line search by default

        // Set up parameters
        LBFGSBParam<double> param = ou.setup_lbfgs_param();

        int number_of_col = global_spline.size();
        int number_of_row = 3; // Should be 3 for XYZ

        // We need to compress g_spline to a single array
        // g_spline -> g_single_array
        // MatrixXd to vector<double>
        // col*0 + (0 to col-1) = x 
        // col*1 + (0 to col-1) = y
        // col*2 + (0 to col-1) = z
        // Initial guess = g_single_array which is x
        Eigen::VectorXd x = 
            VectorXd::Zero(number_of_col * number_of_row);
        x = ou.vect_vector3_to_single(global_spline);

        std::cout << "[optimization_server]" <<
            " x size: " << 
            KBLU << x.size() << KNRM << std::endl;

        // Create solver and function object
        LBFGSBSolver<double> solver(param);  // New solver class
        trajectory_server::optimizer opt(number_of_col * number_of_row);

        // Bounds
        // Setup Lower Bound and upper bound for now we take z to be clamped
        Eigen::VectorXd lb = VectorXd::Zero(number_of_col * number_of_row);
        Eigen::VectorXd ub = VectorXd::Zero(number_of_col * number_of_row);
        for (int i = 0; i < number_of_row; i++)
        {
            for (int j = 0; j < number_of_col; j++)
            {
                if (i != 2)
                {
                    // If representing x and y
                    // Load in max and min double lower and upper bound
                    lb(number_of_col*i + j) = -DBL_MAX;
                    ub(number_of_col*i + j) = DBL_MAX;
                }
                else
                {
                    // If representing z
                    // Load in Lower and Upper bound
                    lb(number_of_col*i + j) = min_height;
                    ub(number_of_col*i + j) = max_height;
                }

            }
        }

        // x will be overwritten to be the best point found
        double fx = 0; // Cost
        opt.set_params_and_data(
            dt, max_acc, protected_zone, reference_spline, obs, time_points);
        opt.set_weights(_weight_vector);
        
        int iter = solver.minimize(opt, x, fx, lb, ub);

        // single_array -- back into --> global_spline
        // vector<double> to vector<Vector3d>

        printf("%s[bspline_optimization.h] Iterations %s%d! %sF(x) %s%lf!\n", 
            KGRN, KNRM, iter, KGRN, KNRM, fx);
        
        vector<Eigen::Vector3d> cp;

        for (int j = 0; j < number_of_col; j++)
        {
            // Add the vector3d to the global_spline vector
            Eigen::Vector3d single_point = Eigen::Vector3d(
                x[number_of_col * 0 + j], x[number_of_col * 1 + j], x[number_of_col * 2 + j]);
            
            cp.push_back(single_point);
        }

        return cp;

    }
}