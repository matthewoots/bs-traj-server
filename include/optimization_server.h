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
#include <random>
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
#include <pcl/kdtree/kdtree_flann.h>

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
            system_clock::time_point origin;
        };
        
        LBFGSBParam<double> setup_lbfgs_param()
        {
            // Set up parameters
            LBFGSBParam<double> param;  // New parameter class
            param.epsilon = 1e-6;
            // param.delta = 1e-3;
            param.max_iterations = 10; // Stable is 20
            param.max_step = 0.20; // default is 1e-4
            param.max_linesearch = 10; // Stable is 20
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

            std::random_device dev;

            vector<trajectory_server::
                optimization_utils::other_agents_traj> other_agents;

            double weight_smooth, weight_feas, weight_term, weight_static, weight_reci;
            
            double fx_smooth, fx_feas, fx_term, fx_static, fx_reci;
            VectorXd grad_smooth, grad_feas, grad_term, grad_static, grad_reci;
            double max_acc, protected_zone, dt;

            vector<Vector3d> reference_cp;
            vector<double> time_points;

            system_clock::time_point origin_time;

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

            void feasibilityCost(vector<Eigen::Vector3d> cp)
            {
                /*
                * Soft limit on the norm of time derivatives 
                * such as velocity, acceleration, jerk and snap.
                */
                double ts_inv2 = 1 / dt / dt;
                double cost = 0;
                int col = (int)cp.size();
                grad_feas = Eigen::VectorXd::Zero(col * 3);

                Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(3,col);
                Eigen::VectorXd grad_single = Eigen::VectorXd::Zero(col * 3);

                /* Check all instances with acceleration limit */
                for (int i = 0; i < cp.size() - 2; i++)
                {
                    Eigen::Vector3d ai = (cp[i + 2] - 2 * cp[i+1] + cp[i]) * ts_inv2;
                    for (int j = 0; j < 3; j++)
                    {
                        if (ai(j) > max_acc)
                        {
                            cost += pow(ai(j) - max_acc, 2);

                            gradient(j, i + 0) += 2 * (ai(j) - max_acc) * ts_inv2;
                            gradient(j, i + 1) += (-4 * (ai(j) - max_acc) * ts_inv2);
                            gradient(j, i + 2) += 2 * (ai(j) - max_acc) * ts_inv2;
                        }
                        else if (ai(j) < -max_acc)
                        {
                            cost += pow(ai(j) + max_acc, 2);

                            gradient(j, i + 0) += 2 * (ai(j) + max_acc) * ts_inv2;
                            gradient(j, i + 1) += (-4 * (ai(j) + max_acc) * ts_inv2);
                            gradient(j, i + 2) += 2 * (ai(j) + max_acc) * ts_inv2;
                        }
                        else {}
                    }
                }
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < col; j++)
                        grad_feas(col * i + j) = gradient(i,j);

                fx_feas = cost;
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
                const double CLEARANCE = protected_zone * 5.0;
                constexpr double a = 1.0, b = 1.0, inv_a2 = 1 / a / a, inv_b2 = 1 / b / b;
                int affected_range = 2;

                for (int i = 0; i < cp.size(); i++)
                {
                    int32_t current_knot_micro = (int32_t)(knots[i] * pow(10,6));
                    system_clock::time_point current_time_on_knot =
                        origin_time + std::chrono::microseconds(current_knot_micro);
                    
                    for (int j = 0; j < other_agents.size(); j++)
                    {
                        if (other_agents[j].id == uav_idx)
                            continue;

                        // double origin_diff = duration<double>(other_agents[j].origin - origin_time).count();
                        // std::cout << "[optimization_server.h] uav " << uav_idx << " on " <<
                        //     other_agents[j].id << " " << KYEL << origin_diff << KNRM << std::endl;
                        
                        std::mt19937 generator(dev());
                        std::uniform_real_distribution<double> dis_middle(-0.3, 0.3);
                        double random_factor = dis_middle(generator);
                        if (random_factor > 0.0 && random_factor < 0.1)
                            random_factor = 0.1;
                        
                        if (random_factor < 0.0 && random_factor > -0.1)
                            random_factor = -0.1;

                        for (int k = 0; k < other_agents[j].cp.size(); k++)
                        {
                            int32_t other_knot_micro = (int32_t)(
                                other_agents[j].knots[k] * pow(10,6));
                            system_clock::time_point other_time_on_knot =
                                other_agents[j].origin + std::chrono::microseconds(
                                other_knot_micro);

                            double time_diff = duration<double>(
                                current_time_on_knot - other_time_on_knot).count();
                            
                            // std::cout << KYEL << "[optimization_server.h] " << time_diff << KNRM << std::endl;

                            if (time_diff > 2.0 * dt)
                            {
                                k += (int)abs(floor(time_diff / dt));
                                continue;
                            }

                            if (time_diff < -2.0 * dt)
                                break;

                            Eigen::Vector3d dist_vec = cp[i] - other_agents[j].cp[k];
                            double smallest_factor = 0.001;
                            if (dist_vec.z() < smallest_factor && dist_vec.z() > -smallest_factor)
                            {
                                dist_vec.z() = dist_vec.z() + random_factor * 0.3;
                            }

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
                                if (i > affected_range && i < cp.size() - affected_range)
                                {
                                    for (int l = 0; l < affected_range*2+1; l++)
                                    {
                                        if (l == affected_range)
                                            continue;
                                        cost = cost + 1 * pow(dist_err, 2);
                                        gradient.col(i - affected_range + l) = 
                                            gradient.col(i - affected_range + l) + 1 * (Coeff.array() * dist_vec.array()).matrix();
                                    }
                                }
                            }
                        }
                        
                    }
                    
                }

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < col; j++)
                        grad_reci(col * i + j) = gradient(i,j);

                fx_reci = cost;

            }

            void staticCollisionCost(vector<Eigen::Vector3d> cp)
            {
                double cost = 0;
                int col = (int)cp.size();
                grad_static = VectorXd::Zero(col * 3);

                MatrixXd gradient = MatrixXd::Zero(3,col);
                VectorXd grad_single = VectorXd::Zero(col * 3);

                double pz_expansion_factor = 2.5;

                for (int i = 0; i < cp.size(); i++)
                { 
                    // For each iteration we should calc this before hand to reduce the load
                    pcl::PointCloud<pcl::PointXYZ>::Ptr tmp_obs(new pcl::PointCloud<pcl::PointXYZ>);
                    tmp_obs = obs_pcl[i];

                    if (tmp_obs->points.size() == 0)
                        continue;

                    vector<Vector3d> kd_points = kdtree_find_points_pcl(cp[i], tmp_obs, 
                        pz_expansion_factor * protected_zone, 10);

                    if (kd_points.size() == 0)
                        continue;

                    // ** Size is too big and optimization may take a few seconds
                    for (int j = 0; j < kd_points.size(); j++)
                    {
                        Vector3d diff = kd_points[j] - cp[i];
                        Vector3d diff_dir = diff / diff.norm();
                        double sq_diff = pow(diff.norm(), 2);
                        double sq_protected = pow(pz_expansion_factor * protected_zone, 2);
                        gradient.col(i) = gradient.col(i) + diff_dir * (sq_protected - sq_diff);
                        cost = cost + (sq_protected - sq_diff);
                        // printf("%s[bspline_optimization.h] cost %lf gradient [%lf %lf %lf]!\n", KBLU, cost, 
                            // gradient.col(i).x(), gradient.col(i).y(), gradient.col(i).z());
                    }               
                    // printf("%s[bspline_optimization.h] updated cost and gradient!\n", KBLU);
                    // printf("%s[bspline_optimization.h] SA : single point Run Time %lf!\n", KCYN, ros::Time::now().toSec() - prev_indi_points);
                }

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < col; j++)
                        grad_static(col * i + j) = gradient(i,j);

                // std::cout << KYEL << "[bspline_optimization.h] S " << grad_static.transpose() << std::endl;

                fx_static = cost;

            }

        public:

            optimizer(int n_) : n(n_) {}

            void set_params_and_data(
                int _uav_idx,
                double _dt, double _max_acc, double _protected_zone,
                vector<Eigen::Vector3d> _reference_cp, 
                pcl::PointCloud<pcl::PointXYZ>::Ptr _full_pcl,
                vector<double> _time_points, 
                vector<trajectory_server::optimization_utils::other_agents_traj> _other_agents,
                system_clock::time_point _origin_time)
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
                origin_time = _origin_time;
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
                feasibilityCost(cp);

                fx = weight_smooth * fx_smooth +
                    weight_feas * fx_feas +
                    // weight_term * fx_term +
                    // weight_static * fx_static +
                    weight_reci * fx_reci;

                for (int i = 0; i < col*3; i++)
                    grad[i] = weight_smooth * grad_smooth(i) +
                        weight_feas * grad_feas(i) + 
                        // weight_term * grad_term(i) +
                        // weight_static * grad_static(i) +
                        weight_reci * grad_reci(i);
                
                return fx;
            }

            vector<Vector3d> kdtree_find_points_pcl(Vector3d point, pcl::PointCloud<pcl::PointXYZ>::Ptr _obs,
            double c, int K)
            {
                pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
                vector<Vector3d> points;

                kdtree.setInputCloud(_obs);

                // Maybe we need to check to see any number here is NA, or inf

                pcl::PointXYZ searchPoint;
                searchPoint.x = point.x();
                searchPoint.y = point.y();
                searchPoint.z = point.z();

                // if (isnan(searchPoint.x) ||
                //     isnan(searchPoint.y) ||
                //     isnan(searchPoint.z) ||
                //     isinf(searchPoint.x) ||
                //     isinf(searchPoint.y) ||
                //     isinf(searchPoint.z) )
                // {
                //     printf("[kdtree_pcl] Found a point with inf or nan!");
                //     std::cout << searchPoint.x << "," << 
                //     searchPoint.y << "," << searchPoint.z << std::endl;
                //     return points;
                // }

                // K nearest neighbor search

                std::vector<int> pointIdxKNNSearch(K);
                std::vector<float> pointKNNSquaredDistance(K);

                // float radius = 256.0f * rand () / (RAND_MAX + 1.0f);

                float radius = (float)c;

                if ( kdtree.nearestKSearch (searchPoint, K, pointIdxKNNSearch, pointKNNSquaredDistance) > 0 )
                {
                    // for (std::size_t i = 0; i < pointIdxRadiusSearch.size (); ++i)
                    for (std::size_t i = 0; i < pointIdxKNNSearch.size (); ++i)
                    {
                        Vector3d kd_point;
                        // When the point is larger than the radius, we do not consider it
                        if (pointKNNSquaredDistance[i] - pow(radius,2) > 0)
                            continue;
                        kd_point.x() = (*_obs)[ pointIdxKNNSearch[i] ].x; 
                        kd_point.y() = (*_obs)[ pointIdxKNNSearch[i] ].y;
                        kd_point.z() = (*_obs)[ pointIdxKNNSearch[i] ].z;
                        points.push_back(kd_point);
                    }
                }

                return points;
            }
    };

};

#endif