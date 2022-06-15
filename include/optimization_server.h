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
using namespace std::chrono;
using namespace std;
using namespace common_utils;
using namespace LBFGSpp;

namespace trajectory_server
{

    class optimization_utils
    {
        public: 

        struct other_agents_traj
        {
            int id;
            vector<Eigen::Vector3d> cp;
            vector<double> knots;
        };
        
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
            int uav_idx;

            vector<trajectory_server::
                optimization_utils::other_agents_traj> other_agents;

            double weight_smooth, weight_feas, weight_term, weight_static, weight_reci;
            
            double fx_smooth, fx_feas, fx_term, fx_static, fx_reci;
            VectorXd grad_smooth, grad_feas, grad_term, grad_static, grad_reci;
            double max_acc, protected_zone, dt;

            vector<Vector3d> reference_cp;
            vector<double> time_points;

            // obs_pcl will save up the points 
            vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> obs_pcl;
            pcl::PointCloud<pcl::PointXYZ>::Ptr full_pcl;

            void smoothnessCost(vector<Eigen::Vector3d> cp)
            {
                /*
                * Minimizing the integral over the squared derivatives 
                * (smoothness) such as acceleration, jerk, snap.
                */
                double cost = 0; 
                int col = (int)cp.size();
                grad_smooth = Eigen::VectorXd::Zero(col * 3);

                // col*0 + (0 to col-1) = x 
                // col*1 + (0 to col-1) = y
                // col*2 + (0 to col-1) = z 

                Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(3,col);
                Eigen::VectorXd grad_single = Eigen::VectorXd::Zero(col * 3);

                /* 3rd derivative, Jerk */
                for (int i = 0; i < (int)cp.size() - 3; i++)
                {
                    Eigen::Vector3d jerk = 
                        cp[i + 3] - 3 * cp[i + 2] + 3 * cp[i + 1] - cp[i];

                    cost = cost + pow(jerk.norm(), 2);
                    Eigen::Vector3d tmp_j = 2.0 * jerk;

                    gradient.col(i + 0) += (-tmp_j);
                    gradient.col(i + 1) += (3.0 * tmp_j);
                    gradient.col(i + 2) += (-3.0 * tmp_j);
                    gradient.col(i + 3) += (tmp_j);
                }
                
                /* 2nd derivative, Acceleration */
                for (int i = 0; i < (int)cp.size() - 2; i++)
                {   
                    Eigen::Vector3d acc = 
                        cp[i + 2] - 2 * cp[i + 1] + cp[i];

                    cost = cost + pow(acc.norm(), 2);
                    Vector3d tmp_a = 2.0 * acc;
                    
                    gradient.col(i + 0) = (tmp_a);
                    gradient.col(i + 1) = (-2.0 * tmp_a);
                    gradient.col(i + 2) = (tmp_a);
                }

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < col; j++)
                        grad_smooth(col * i + j) = gradient(i,j);
                
                fx_smooth = cost;
            }

            void reciprocalAvoidanceCost(vector<Eigen::Vector3d> cp, vector<double> knots)
            {
                /*
                * Reciprocal avoidance function that handles 
                * a collision-free path between agents
                */
                double cost = 0;
                int col = (int)cp.size();
                double tolerance = knots[1] - knots[0];
                // std::cout << KGRN << "[bspline_optimization.h] tolerance " << tolerance << std::endl;

                grad_reci = VectorXd::Zero(col * 3);

                MatrixXd gradient = MatrixXd::Zero(3,col);
                VectorXd grad_single = VectorXd::Zero(col * 3);

                double magnitude = 1;
                const double CLEARANCE = protected_zone * 2.0;
                constexpr double a = 2.0, b = 1.0, inv_a2 = 1 / a / a, inv_b2 = 1 / b / b;

                for (int i = 0; i < cp.size(); i++)
                {
                    // std::cout << KYEL << "[bspline_optimization.h] knots " << i << " " << knots(i) << std::endl;
                    double current_knot = knots[i];
                    for (int j = 0; j < other_agents.size(); j++)
                    {
                        if (other_agents[j].id == uav_idx)
                            continue;
                        for (int k = 0; k < other_agents[j].cp.size(); k++)
                        {
                            Eigen::Vector3d dist_vec = cp[i] - other_agents[j].cp[k];
                            double ellip_dist = sqrt(dist_vec(2) * dist_vec(2) * inv_a2 + (dist_vec(0) * dist_vec(0) + dist_vec(1) * dist_vec(1)) * inv_b2);
                            
                            double dist_err = CLEARANCE - ellip_dist;

                            Vector3d Coeff = Vector3d(
                                -2 * (CLEARANCE / ellip_dist - 1) * inv_b2,
                                -2 * (CLEARANCE / ellip_dist - 1) * inv_b2,
                                -2 * (CLEARANCE / ellip_dist - 1) * inv_a2
                            );

                            if (dist_err < 0 || dist_vec.norm() == 0)
                            {
                                // Empty
                            }
                            else
                            {
                                // std::cout << KYEL << "[bspline_optimization.h] dist_err " << dist_err 
                                //     << " magnitude " << magnitude << " dist_vec " << dist_vec.norm() << std::endl;
                                cost = cost + magnitude * pow(dist_err, 2);
                                gradient.col(i) = gradient.col(i) + magnitude * (Coeff.array() * dist_vec.array()).matrix();
                            }
                        }
                        
                    }
                    
                }

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < col; j++)
                        grad_reci(col * i + j) = gradient(i,j);

                fx_reci = cost;

            }

        public:

            optimizer(int n_) : n(n_) {}

            void set_params_and_data(
                int _uav_idx,
                double _dt, double _max_acc, double _protected_zone,
                vector<Eigen::Vector3d> _reference_cp, 
                pcl::PointCloud<pcl::PointXYZ>::Ptr _full_pcl,
                vector<double> _time_points, 
                vector<trajectory_server::optimization_utils::other_agents_traj> _other_agents)
            {
                uav_idx = _uav_idx;
                dt = _dt; 
                max_acc = _max_acc;
                protected_zone = _protected_zone;

                // Reset values
                reference_cp.clear();
                _full_pcl->points.clear();

                reference_cp = _reference_cp;
                full_pcl = _full_pcl;
                time_points = _time_points;
                other_agents = _other_agents;
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
                double factor = 2.0;
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

                smoothnessCost(cp);
                reciprocalAvoidanceCost(cp, time_points);

                fx = weight_smooth * fx_smooth +
                    // weight_feas * fx_feas +
                    // weight_term * fx_term +
                    // weight_static * fx_static +
                    weight_reci * fx_reci;

                for (int i = 0; i < col*3; i++)
                    grad[i] = weight_smooth * grad_smooth(i) +
                        // weight_feas * grad_feas(i) + 
                        // weight_term * grad_term(i) +
                        // weight_static * grad_static(i) +
                        weight_reci * grad_reci(i);
                
                return fx;
            }
    };

};

#endif