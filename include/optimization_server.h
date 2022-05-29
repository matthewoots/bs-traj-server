/*
* optimization_server.h
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

#ifndef OPTIMIZATION_SERVER_H
#define OPTIMIZATION_SERVER_H

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

#include "bspline_utils.hpp"
#include "trajectory_server.h"

#include <pcl/point_types.h>
#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/crop_box.h>

#include "LBFGSB.h"
#include "common_utils.h"

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
using namespace common_utils;
using namespace LBFGSpp;

namespace trajectory_server
{

    class optimization_utils
    {
        public: 
        
        LBFGSBParam<double> setup_lbfgs_param()
        {
            // Set up parameters
            LBFGSBParam<double> param;  // New parameter class
            param.epsilon = 1e-6;
            // param.delta = 1e-3;
            param.max_iterations = 20; // Stable is 20
            param.max_step = 0.20; // default is 1e-4
            param.max_linesearch = 20; // Stable is 20
            // param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
        
            return param;
        }

        Eigen::VectorXd vect_vector3_to_single(vector<Eigen::Vector3d> vector3)
        {
            // We need to compress g_spline to a single array
            // g_spline -> g_single_array
            // MatrixXd to vector<double>
            // col*0 + (0 to col-1) = x 
            // col*1 + (0 to col-1) = y
            // col*2 + (0 to col-1) = z
            int number_of_col = vector3.size();
            Eigen::VectorXd single_array = 
                VectorXd::Zero(number_of_col * 3);
            
            for (int j = 0; j < number_of_col; j++)
            {
                single_array(number_of_col * 0 + j) = vector3[j][0];
                single_array(number_of_col * 1 + j) = vector3[j][1];
                single_array(number_of_col * 2 + j) = vector3[j][2];
            }
        
            return single_array;
        }
        
    };

    class optimizer
    {
        private:

            int n;

            double weight_smooth, weight_feas, weight_term, weight_static, weight_reci;
            double max_acc, protected_zone, dt;

            vector<Vector3d> reference_cp;
            vector<double> time_points;

            // obs_pcl will save up the points 
            vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> obs_pcl;
            pcl::PointCloud<pcl::PointXYZ>::Ptr full_pcl;

        public:
            optimizer(int n_) : n(n_) {}

            void set_params_and_data(
                double _dt, double _max_acc, double _protected_zone,
                vector<Eigen::Vector3d> _reference_cp, 
                pcl::PointCloud<pcl::PointXYZ>::Ptr _full_pcl,
                vector<double> _time_points)
            {
                dt = _dt; 
                max_acc = _max_acc;
                protected_zone = _protected_zone;

                // Reset values
                reference_cp.clear();
                _full_pcl->points.clear();

                reference_cp = _reference_cp;
                full_pcl = _full_pcl;
                time_points = _time_points;
            }

            void set_weights(vector<double> weight_vector)
            {
                // Index for weights
                // 0. _weight_smooth = weight_vector[0]; 
                // 1. _weight_feas = weight_vector[1];
                // 2. _weight_term = weight_vector[2];
                // 3. _weight_static = weight_vector[3];
                // 4. _weight_reci = weight_vector[4];
                weight_smooth = weight_vector[0]; 
                weight_feas = weight_vector[1];
                weight_term = weight_vector[2];
                weight_static = weight_vector[3];
                weight_reci = weight_vector[4];
            }

            double operator()(const VectorXd& x, VectorXd& grad)
            {
                double fx = 0.0;            
                obs_pcl.clear();
                double factor = 4.0;
                Eigen::Vector3d dimensions = Eigen::Vector3d(
                    factor * protected_zone,
                    factor * protected_zone,
                    factor * protected_zone);

                int col = (int)x.size() / 3;

                // single_array -- into --> global_spline
                // vector<double> to vector<Vector3d>
                vector<Eigen::Vector3d> cp;
                for (int j = 0; j < col; j++)
                {
                    // Add the vector3d to the global_spline vector
                    Eigen::Vector3d single_point = Eigen::Vector3d(
                        x[col * 0 + j], x[col * 1 + j], x[col * 2 + j]);
                    
                    cp.push_back(single_point);

                    // Add the localized pcl to the obstacle vector
                    pcl::PointCloud<pcl::PointXYZ>::Ptr _obs(
                        new pcl::PointCloud<pcl::PointXYZ>);
                    
                    _obs = common_utils::pcl_ptr_box_crop(
                        full_pcl, single_point, dimensions);

                    obs_pcl.push_back(_obs);
                }
                
                return fx;
            }
    };

    class optimization_server
    {
        public:
        
            optimization_server(vector<double> weight_vector)
            {
                if (weight_vector.size() != 5)
                    return;
                // Index for weights
                // 0. _weight_smooth = weight_vector[0]; 
                // 1. _weight_feas = weight_vector[1];
                // 2. _weight_term = weight_vector[2];
                // 3. _weight_static = weight_vector[3];
                // 4. _weight_reci = weight_vector[4];
                _weight_vector = weight_vector;
            }
            
            ~optimization_server(){}

            // Solver will solve with the given partial global_spline 
            vector<Eigen::Vector3d> solver(
                vector<Eigen::Vector3d> global_spline, 
                vector<Eigen::Vector3d> reference_spline, 
                vector<double> time_points, 
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs, 
                double dt, double max_acc,
                double protected_zone,
                double min_height, double max_height);


        private:

            trajectory_server::optimization_utils ou;

            vector<double> _weight_vector;

    };

};

#endif