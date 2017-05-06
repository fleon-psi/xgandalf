/*
 * pointAutocorrelation.h
 *
 *  Created on: 06.05.2017
 *      Author: Yaro
 */

#ifndef POINTAUTOCORRELATION_H_
#define POINTAUTOCORRELATION_H_

#include <Eigen/Dense> 
#include <limits>

// all autocorrelation results will have only half of the possible points, since symmetric points (at x < 0) will be removed

void getPointAutocorrelation(Eigen::Matrix3Xf& autocorrelationPoints, const Eigen::Matrix3Xf& points,
        float minNormInAutocorrelation, float maxNormInAutocorrelation);

void getPointAutocorrelation(Eigen::Matrix3Xf& autocorrelationPoints, Eigen::VectorXf& centerPointIndices, Eigen::VectorXf& shiftedPointIndices,
        const Eigen::Matrix3Xf& points, float minNormInAutocorrelation, float maxNormInAutocorrelation);

#endif /* POINTAUTOCORRELATION_H_ */
