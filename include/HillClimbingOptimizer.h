/*
 * HillClimbingOptimizer.h
 *
 *  Created on: 16.04.2017
 *      Author: Yaro
 */

#ifndef HILLCLIMBINGOPTIMIZER_H_
#define HILLCLIMBINGOPTIMIZER_H_

#include <InverseSpaceTransform.h>
#include <Eigen/Dense>

class HillClimbingOptimizer {
public:
    HillClimbingOptimizer();

    //watch out! gradient, closeToPeaksCount and inverseTransformEvaluation are changed in this function (for performance reasons)!
    void computeStep(Eigen::Matrix3Xf& gradient, Eigen::RowVectorXf& closeToPeaksCount, Eigen::RowVectorXf& inverseTransformEvaluation,
            bool useStepOrthogonalizationFlag);

private:
    //accuracy constants
    InverseSpaceTransform *transform;

    float gamma;
    int initialIterationCount;
    float maxStep;
    float minStep;
    int calmDownIterationCount;
    float calmDownFactor;
    float directionChangeFactor;
    int localFitIterationCount;
    int localCalmDownIterationCount;
    float localCalmDownFactor;

    // interna
    Eigen::Matrix3Xf stepDirection;
    Eigen::Array< float, 1, Eigen::Dynamic > stepLength;
    Eigen::Matrix3Xf step;

    Eigen::Matrix3Xf previousStepDirection;
    Eigen::Array< float, 1, Eigen::Dynamic > previousStepLength;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* HILLCLIMBINGOPTIMIZER_H_ */
