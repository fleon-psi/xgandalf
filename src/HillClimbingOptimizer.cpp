/*
 * HillClimbingOptimizer.cpp
 *
 *  Created on: 16.04.2017
 *      Author: Yaro
 */

#include <HillClimbingOptimizer.h>
#include <cstddef>

using namespace Eigen;
using namespace std;

HillClimbingOptimizer::HillClimbingOptimizer() :
        transform(NULL), gamma(0), initialIterationCount(0), maxStep(0), minStep(0), calmDownIterationCount(0),
                calmDownFactor(0), directionChangeFactor(0), localFitIterationCount(0),
                localCalmDownIterationCount(0), localCalmDownFactor(0)
{
}

void HillClimbingOptimizer::computeStep(Matrix3Xf& gradient, RowVectorXf& closeToPeaksCount, RowVectorXf& inverseTransformEvaluation,
        bool useStepOrthogonalization)
{
    //reuse memory for processing for performance reasons
    Matrix3Xf& stepDirection = gradient;
    RowVectorXf& closeToPeaksFactor = closeToPeaksCount;
    RowVectorXf& functionEvaluationFactor = inverseTransformEvaluation;

    stepDirection = gradient.colwise().normalized();

    Array< float, 1, Dynamic > directionChange = (stepDirection.array() * previousStepDirection.array()).matrix().colwise().squaredNorm();

    Array< float, 1, Dynamic > stepDirectionFactor = ((directionChange + 1) / 2).square().square() * 1.5 + 0.5;   //directionChange in [-1 1]   
    closeToPeaksFactor = (-1 * closeToPeaksCount.array() + 0.8).square() * 6 + 0.5;   //closeToPeaks in [0 1]
    functionEvaluationFactor = ((-1 * inverseTransformEvaluation.array() + 0.8) / 2).cube() * 4 + 0.3;   //functionEvaluation in [-1 1]

    if (useStepOrthogonalization) {
        for (int i = 0; i < closeToPeaksCount.size(); i++) {
            if (directionChange(i) < -0.4) {
                stepDirection.col(i) = (stepDirection.col(i) + previousStepDirection.col(i)).normalized();
                stepDirectionFactor(i) = directionChangeFactor;
            }
        }
    }

    stepLength = ((0.5 * (minStep + (maxStep - minStep) * gamma) + 0.5 * previousStepLength.array())
            * stepDirectionFactor * closeToPeaksFactor.array() * functionEvaluationFactor.array()).min(maxStep).max(minStep);
    
    step = stepDirection.array().rowwise()*stepLength;
}
