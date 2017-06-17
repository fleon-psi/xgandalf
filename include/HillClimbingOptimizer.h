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
    typedef struct {
        float gamma;
        float maxStep;
        float minStep;
        float directionChangeFactor;
    } stepComputationAccuracyConstants_t;

    typedef struct {
        int initialIterationCount;
        int calmDownIterationCount;
        float calmDownFactor;
        int localFitIterationCount;
        int localCalmDownIterationCount;
        float localCalmDownFactor;

        int functionSelection;
        float optionalFunctionArgument;
        float maxCloseToPointDeviation;

        stepComputationAccuracyConstants_t stepComputationAccuracyConstants;
    } hillClimbingAccuracyConstants_t;

    HillClimbingOptimizer();

    void performOptimization(const Eigen::Matrix3Xf& pointsToTransform, Eigen::Matrix3Xf& positionsToOptimize);
    Eigen::RowVectorXf& getLastInverseTransformEvaluation();
    Eigen::RowVectorXf& getCloseToPointsCount();
    std::vector< std::vector< uint16_t > >& getPointsCloseToEvaluationPositions_indices();

    void setHillClimbingAccuracyConstants(hillClimbingAccuracyConstants_t accuracyConstants);

    //optional
    void setPointsToTransformWeights(const Eigen::RowVectorXf& pointsToTransformWeights);

public:
    void setStepComputationAccuracyConstants(stepComputationAccuracyConstants_t stepComputationAccuracyConstants);

    //watch out! gradient, closeToPointsCount and inverseTransformEvaluation are changed in this function (for performance reasons)!
    void computeStep(Eigen::Matrix3Xf& gradient, Eigen::RowVectorXf& closeToPointsCount, Eigen::RowVectorXf& inverseTransformEvaluation,
            bool useStepOrthogonalization);

    void performOptimizationStep(Eigen::Matrix3Xf& positionsToOptimize, bool useStepOrthogonalization);

    InverseSpaceTransform transform;
    hillClimbingAccuracyConstants_t hillClimbingAccuracyConstants;

    // interna
    Eigen::Matrix3Xf step;

    Eigen::Matrix3Xf previousStepDirection;
    Eigen::Array< float, 1, Eigen::Dynamic > previousStepLength;

    Eigen::RowVectorXf lastInverseTransformEvaluation;
};

#endif /* HILLCLIMBINGOPTIMIZER_H_ */
