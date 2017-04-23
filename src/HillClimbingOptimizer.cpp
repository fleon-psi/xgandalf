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
        transform(), initialIterationCount(0), calmDownIterationCount(0), calmDownFactor(0), localFitIterationCount(0),
                localCalmDownIterationCount(0), localCalmDownFactor(0), gamma(0), maxStep(0), minStep(0), directionChangeFactor(0)
{
}

void HillClimbingOptimizer::performOptimization(const Matrix3Xf& pointsToTransform, Matrix3Xf& positionsToOptimize)
{
//    std::ofstream ofs("workfolder/tmp", std::ofstream::out);
//    ofs << positionsToOptimize.transpose().eval() << endl;

    transform.setPointsToTransform(pointsToTransform);

    float gamma_initial = this->gamma;
    float maxStep_initial = this->maxStep;
    float minStep_initial = this->minStep;

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

    for (int i = 0; i < initialIterationCount; i++) {
        transform.setLocalTransformFlag();
        transform.clearRadialWeightingFlag();
        bool useStepOrthogonalization = true;

        performOptimizationStep(positionsToOptimize, useStepOrthogonalization);
//        ofs << positionsToOptimize.transpose().eval() << endl;
    }

    for (int i = 0; i < calmDownIterationCount; i++) {
        transform.setLocalTransformFlag();
        transform.clearRadialWeightingFlag();
        bool useStepOrthogonalization = false;

        maxStep = maxStep * localCalmDownFactor;
        minStep = minStep * localCalmDownFactor;
        gamma = gamma * localCalmDownFactor;

        performOptimizationStep(positionsToOptimize, useStepOrthogonalization);
//        ofs << positionsToOptimize.transpose().eval() << endl;
    }

    this->gamma = gamma_initial;
    this->maxStep = maxStep_initial;
    this->minStep = minStep_initial;

    // can be optimized! Does not always need to compute slope, closeToPeaks and gradient
    transform.performTransform(positionsToOptimize);
    lastInverseTransformEvaluation = transform.getInverseTransformEvaluation();
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

const RowVectorXf& HillClimbingOptimizer::getLastInverseTransformEvaluation()
{
    return lastInverseTransformEvaluation;
}

const RowVectorXf& HillClimbingOptimizer::getCloseToPeaksCount()
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
                stepDirectionFactor(i) = directionChangeFactor;
            }
        }
    }
//    cout << "stepDirection " << endl << stepDirection << endl << endl;

    previousStepDirection = stepDirection;

    Array< float, 1, Eigen::Dynamic >& stepLength = previousStepLength;  //reuse memory
    stepLength = ((0.5 * (minStep + (maxStep - minStep) * gamma) + 0.5 * previousStepLength)
            * stepDirectionFactor * closeToPeaksFactor.array() * functionEvaluationFactor.array()).min(maxStep).max(minStep);

    step = stepDirection.array().rowwise() * stepLength;
}

void HillClimbingOptimizer::setHillClimbingStrategyAccuracyConstants(int initialIterationCount, int calmDownIterationCount, float calmDownFactor,
        int localFitIterationCount, int localCalmDownIterationCount, float localCalmDownFactor)
{
    this->initialIterationCount = initialIterationCount;
    this->calmDownIterationCount = calmDownIterationCount;
    this->calmDownFactor = calmDownFactor;
    this->localFitIterationCount = localFitIterationCount;
    this->localCalmDownIterationCount = localCalmDownIterationCount;
    this->localCalmDownFactor = localCalmDownFactor;

}
void HillClimbingOptimizer::setStepComputationAccuracyConstants(float gamma, float minStep, float maxStep, float directionChangeFactor)
{
    this->gamma = gamma;
    this->minStep = minStep;
    this->maxStep = maxStep;
    this->directionChangeFactor = directionChangeFactor;
}

void HillClimbingOptimizer::setInverseSpaceTransformAccuracyConstants(int functionSelection, float optionalFunctionArgument, float maxCloseToPeakDeviation)
{
    transform.setFunctionSelection(functionSelection);
    transform.setOptionalFunctionArgument(optionalFunctionArgument);
    transform.setMaxCloseToPeakDeviation(maxCloseToPeakDeviation);
}
