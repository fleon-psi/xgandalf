/*
 * HillClimbingOptimizer.cpp
 *
 *  Created on: 16.04.2017
 *      Author: Yaro
 */

#include <HillClimbingOptimizer.h>
#include <cstddef>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

HillClimbingOptimizer::HillClimbingOptimizer() :
        transform(), hillClimbingAccuracyConstants()
{
}

void HillClimbingOptimizer::performOptimization(const Matrix3Xf& pointsToTransform, Matrix3Xf& positionsToOptimize)
{
//    std::ofstream ofs("workfolder/tmp", std::ofstream::out);
//    ofs << positionsToOptimize.transpose().eval() << endl;

    transform.setPointsToTransform(pointsToTransform);

    float& gamma = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.gamma;
    float& maxStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.maxStep;
    float& minStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.minStep;

    float gamma_initial = gamma;
    float maxStep_initial = maxStep;
    float minStep_initial = minStep;

    const int initialIterationCount = hillClimbingAccuracyConstants.initialIterationCount;
    const int calmDownIterationCount = hillClimbingAccuracyConstants.calmDownIterationCount;
    const float calmDownFactor = hillClimbingAccuracyConstants.calmDownFactor;
    const int localFitIterationCount = hillClimbingAccuracyConstants.localFitIterationCount;
    const int localCalmDownIterationCount = hillClimbingAccuracyConstants.localCalmDownIterationCount;
    const float localCalmDownFactor = hillClimbingAccuracyConstants.localCalmDownFactor;

    previousStepDirection = Matrix3Xf::Zero(3, positionsToOptimize.cols());
    previousStepLength = Array< float, 1, Eigen::Dynamic >::Constant(1, positionsToOptimize.cols(), minStep + (maxStep - minStep) / 4);

    for (int i = 0; i < initialIterationCount; i++) {
        transform.clearLocalTransformFlag();
        transform.setRadialWeightingFlag();
        bool useStepOrthogonalization = true;

        performOptimizationStep(positionsToOptimize, useStepOrthogonalization);
//        ofs << positionsToOptimize.transpose().eval() << endl;
    }

    for (int i = 0; i < calmDownIterationCount; i++) {
        transform.clearLocalTransformFlag();
        transform.setRadialWeightingFlag();
        bool useStepOrthogonalization = true;

        maxStep = maxStep * calmDownFactor;
        minStep = minStep * calmDownFactor;
        gamma = gamma * calmDownFactor;

        performOptimizationStep(positionsToOptimize, useStepOrthogonalization);
//        ofs << positionsToOptimize.transpose().eval() << endl;
    }

    for (int i = 0; i < localFitIterationCount; i++) {
        transform.setLocalTransformFlag();
        transform.clearRadialWeightingFlag();
        bool useStepOrthogonalization = true;

        performOptimizationStep(positionsToOptimize, useStepOrthogonalization);
//        ofs << positionsToOptimize.transpose().eval() << endl;
    }

    for (int i = 0; i < localCalmDownIterationCount; i++) {
        transform.setLocalTransformFlag();
        transform.clearRadialWeightingFlag();
        bool useStepOrthogonalization = false;

        maxStep = maxStep * localCalmDownFactor;
        minStep = minStep * localCalmDownFactor;
        gamma = gamma * localCalmDownFactor;

        performOptimizationStep(positionsToOptimize, useStepOrthogonalization);
//        ofs << positionsToOptimize.transpose().eval() << endl;
    }

    // can be optimized! Does not always need to compute slope, closeToPeaks and gradient
    transform.performTransform(positionsToOptimize);
    lastInverseTransformEvaluation = transform.getInverseTransformEvaluation();

    gamma = gamma_initial;
    maxStep = maxStep_initial;
    minStep = minStep_initial;
}

void HillClimbingOptimizer::performOptimizationStep(Matrix3Xf& positionsToOptimize, bool useStepOrthogonalization)
{
    transform.performTransform(positionsToOptimize);
    Matrix3Xf& gradient = transform.getGradient();
    RowVectorXf& closeToPeaksCount = transform.getCloseToPointsCount();
    RowVectorXf& inverseTransformEvaluation = transform.getInverseTransformEvaluation();
    computeStep(gradient, closeToPeaksCount, inverseTransformEvaluation, useStepOrthogonalization);
//    cout << "step: " << endl << step << endl << endl;
    positionsToOptimize += step;
//    cout << "positionsToOptimize: " << endl << positionsToOptimize << endl << endl;
}

void HillClimbingOptimizer::setPointsToTransformWeights(const Eigen::RowVectorXf pointsToTransformWeights)
{
    transform.setPointsToTransformWeights(pointsToTransformWeights);
}

RowVectorXf& HillClimbingOptimizer::getLastInverseTransformEvaluation()
{
    return lastInverseTransformEvaluation;
}

RowVectorXf& HillClimbingOptimizer::getCloseToPeaksCount()
{
    return transform.getCloseToPointsCount();
}

vector< vector< uint16_t > >& HillClimbingOptimizer::getPeaksCloseToEvaluationPositions_indices()
{
    return transform.getPointsCloseToEvaluationPositions_indices();
}

void HillClimbingOptimizer::computeStep(Matrix3Xf& gradient, RowVectorXf& closeToPeaksCount, RowVectorXf& inverseTransformEvaluation,
        bool useStepOrthogonalization)
{
//reuse memory for processing for performance reasons
    Matrix3Xf& stepDirection = gradient;
    RowVectorXf& closeToPeaksFactor = closeToPeaksCount;
    RowVectorXf& functionEvaluationFactor = inverseTransformEvaluation;

    stepDirection = gradient.colwise().normalized();
//    cout << "stepDirection " << endl << stepDirection << endl << endl<< previousStepDirection << endl << endl;

    Array< float, 1, Dynamic > directionChange = (stepDirection.array() * previousStepDirection.array()).matrix().colwise().sum();
//    cout << "dirChange " << endl << directionChange << endl << endl;

    Array< float, 1, Dynamic > stepDirectionFactor = ((directionChange + 1) / 2).square().square() * 1.5 + 0.5;   //directionChange in [-1 1]   
    closeToPeaksFactor = (-1 * closeToPeaksCount.array() + 0.8).square() * 6 + 0.5;   //closeToPeaks in [0 1]
    functionEvaluationFactor = ((-1 * inverseTransformEvaluation.array() + 0.8) / 2).cube() * 4 + 0.3;   //functionEvaluation in [-1 1]
//    cout << stepDirectionFactor << endl << endl << closeToPeaksFactor << endl << endl << functionEvaluationFactor << endl;

//    cout << "stepDirection " << endl << stepDirection << endl << endl;
    if (useStepOrthogonalization) {
        for (int i = 0; i < closeToPeaksCount.size(); i++) {
            if (directionChange(i) < -0.4) {
                stepDirection.col(i) = (stepDirection.col(i) + previousStepDirection.col(i)).normalized();
                stepDirectionFactor(i) = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.directionChangeFactor;
            }
        }
    }
//    cout << "stepDirection " << endl << stepDirection << endl << endl;

    previousStepDirection = stepDirection;

    const float minStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.minStep;
    const float maxStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.maxStep;
    const float gamma = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.gamma;

    Array< float, 1, Eigen::Dynamic >& stepLength = previousStepLength;  //reuse memory
    stepLength = ((0.5 * (minStep + (maxStep - minStep) * gamma) + 0.5 * previousStepLength)
            * stepDirectionFactor * closeToPeaksFactor.array() * functionEvaluationFactor.array()).min(maxStep).max(minStep);

    step = stepDirection.array().rowwise() * stepLength;
}

void HillClimbingOptimizer::setStepComputationAccuracyConstants(stepComputationAccuracyConstants_t stepComputationAccuracyConstants)
{
    hillClimbingAccuracyConstants.stepComputationAccuracyConstants = stepComputationAccuracyConstants;
}

void HillClimbingOptimizer::setHillClimbingAccuracyConstants(hillClimbingAccuracyConstants_t accuracyConstants)
{
    hillClimbingAccuracyConstants = accuracyConstants;

    transform.setFunctionSelection(accuracyConstants.functionSelection);
    transform.setOptionalFunctionArgument(accuracyConstants.optionalFunctionArgument);
    transform.setMaxCloseToPeakDeviation(accuracyConstants.maxCloseToPeakDeviation);
}
