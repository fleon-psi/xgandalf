/*
 * eigenSTLContainers.h
 *
 *  Created on: 21.12.2016
 *      Author: Yaro
 */

#ifndef EIGENSTLCONTAINERS_H_
#define EIGENSTLCONTAINERS_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <vector>

namespace EigenSTL
{
    typedef std::vector<Eigen::Vector2f, Eigen::aligned_allocator<Eigen::Vector2f>> vector_Vector2f;
    typedef std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> vector_Vector2d;
    typedef std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> vector_Vector3f;
    typedef std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> vector_Vector3i;
    typedef std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vector_Vector3d;
    typedef std::vector<Eigen::Vector4f, Eigen::aligned_allocator<Eigen::Vector4f>> vector_Vector4f;

    typedef std::vector<Eigen::ArrayXXd, Eigen::aligned_allocator<Eigen::ArrayXXd>> vector_ArrayXXd;

    typedef std::vector<Eigen::Matrix3Xd, Eigen::aligned_allocator<Eigen::Matrix3Xd>> vector_Matrix3Xd;

    typedef std::vector<Eigen::RowVectorXd, Eigen::aligned_allocator<Eigen::RowVectorXd>> vector_RowVectorXd;
}

#endif /* EIGENSTLCONTAINERS_H_ */
