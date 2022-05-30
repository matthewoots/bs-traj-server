/*
 * common_utils.h
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
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 * 
 * 
 */

#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <Eigen/Dense>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/crop_box.h>

using namespace Eigen;
using namespace std;

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

namespace common_utils
{
    /** @brief Filter/Crop point cloud with the dimensions given */
    inline pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_ptr_box_crop(
        pcl::PointCloud<pcl::PointXYZ>::Ptr _pc, 
        Eigen::Vector3d centroid, Eigen::Vector3d dimension)
    {   
        pcl::PointCloud<pcl::PointXYZ>::Ptr output(
            new pcl::PointCloud<pcl::PointXYZ>);

        Eigen::Vector3d min = centroid - (dimension / 2);
        Eigen::Vector3d max = centroid + (dimension / 2);

        pcl::CropBox<pcl::PointXYZ> box_filter;
        box_filter.setMin(Eigen::Vector4f(min.x(), min.y(), min.z(), 1.0));
        box_filter.setMax(Eigen::Vector4f(max.x(), max.y(), max.z(), 1.0));

        box_filter.setInputCloud(_pc);
        box_filter.filter(*output);

        return output;
    }
};

#endif