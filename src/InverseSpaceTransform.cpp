/*
 * InverseSpaceTransform.cpp
 *
 *  Created on: 11.04.2017
 *      Author: Yaro
 */

#include <InverseSpaceTransform.h>
#include <assert.h>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace Eigen;
using namespace std;

static inline void function1(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument);
static inline void function2(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
static inline void function3(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
static inline void function4(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
static inline void function5(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
static inline void function6(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
static inline void function7(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument);
static inline void function8(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
static inline void function9(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument);

InverseSpaceTransform::InverseSpaceTransform() :
        functionSelection(0), optionalFunctionArgument(1), localTransform(false), radialWeighting(false), maxCloseToPeakDeviation(0.15), resultsUpToDate(false)
{

}

InverseSpaceTransform::InverseSpaceTransform(float maxCloseToPeakDeviation) :
        functionSelection(0), optionalFunctionArgument(1), localTransform(false), radialWeighting(false), maxCloseToPeakDeviation(maxCloseToPeakDeviation),
                resultsUpToDate(false)
{

}

void InverseSpaceTransform::performTransform(const Matrix3Xf& positionsToEvaluate)
{
    float pointsToTransformCount_inverse = 1 / (float) pointsToTransform.cols();

    ArrayXXf x = pointsToTransform.transpose() * positionsToEvaluate;
    onePeriodicFunction(x);

//    cout << slope << endl << endl << pointsToTransform << endl << endl << pointsToTransformWeights << endl << endl;
    Matrix3Xf fullGradient;
    if (localTransform) {
        fullGradient = (pointsToTransform.array().rowwise() * pointsToTransformWeights.array()).matrix() * slope.matrix();
        functionEvaluation = functionEvaluation * closeToPeak.matrix().cast< float >().array();
        slope = slope * closeToPeak.matrix().cast< float >().array();
    }
//    cout << slope << endl << endl << functionEvaluation << endl << endl << fullGradient << endl << endl;

    gradient = (pointsToTransform.array().rowwise() * pointsToTransformWeights.array()).matrix() * slope.matrix();
    inverseTransformEvaluation = pointsToTransformWeights * functionEvaluation.matrix() * (1 / inverseTransformEvaluationScalingFactor);

    if (closeToPeak.rows() <= 255) {
        closeToPeaksCount = closeToPeak.matrix().cast< uint8_t >().colwise().sum().cast< float >();
    } else {
        closeToPeaksCount = closeToPeak.matrix().cast< uint16_t >().colwise().sum().cast< float >();
    }
//    cout << gradient << endl << endl << inverseTransformEvaluation << endl << endl << closeToPeak << endl << endl;

    if (localTransform) {
//        gradient = gradient.array().rowwise() / closeToPeaksCount.array();
//        gradient = (gradient.array() == 0).select(fullGradient * pointsToTransformCount_inverse, gradient);
        for (int i = 0; i < closeToPeaksCount.size(); i++) {
            gradient.col(i) =
                    (closeToPeaksCount(i) != 0) ? (gradient.col(i) * (1.0f / closeToPeaksCount(i))) : (fullGradient.col(i) * pointsToTransformCount_inverse);
        }
    } else {
        gradient = gradient * pointsToTransformCount_inverse;
    }

    closeToPeaksCount = closeToPeaksCount * pointsToTransformCount_inverse;

    resultsUpToDate = true;
}

// cos(x * 2*pi).^optionalFunctionArgument
static inline void function1(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument)
{
    assert((optionalFunctionArgument - round(optionalFunctionArgument)) == 0);
    assert((int )optionalFunctionArgument % 2 != 0); // An even optionanFunctionArgument does not make sense with this function.

    if (optionalFunctionArgument == 1) {
        functionEvaluation = cos(2 * M_PI * x);
        slope = -sin(2 * M_PI * x);
    } else {
        functionEvaluation = pow(cos(2 * M_PI * x), (int) optionalFunctionArgument);
        float n = optionalFunctionArgument;
        float scaling = pow(n, n / 2) / pow((n - 1), (int) (n - 1) / 2);
        slope = -scaling * sin(2 * M_PI * x) * pow(cos(2 * M_PI * x), (int) (n - 1));
    }
}

// sawtooth
static inline void function2(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
{
    functionEvaluation = (-4 * abs(x) + 1);
    slope = -x.sign();
}

static inline void function3(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
{
    auto x_p4 = x.square().square();

    functionEvaluation = (32.4645 * x_p4 - 16.2782 * x.square() + 1);
    slope = 20.711 * x.cube() - 5.19241 * x;
}

//amplitude from 2, slope from 3
static inline void function4(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
{
    functionEvaluation = (-4 * abs(x) + 1);
    slope = 20.711 * x.cube() - 5.19241 * x;
}

static inline void function5(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
{
    functionEvaluation = -8 * x.square() + 1;
    slope = -2 * x;
}

static inline void function6(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
{
    functionEvaluation = -32 * x.square().square() + 1;
    slope = -8 * x.cube();
}

static inline void function7(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument)
{
    functionEvaluation = -(abs(x) - optionalFunctionArgument / 2).sign();
    slope = -x.sign() * (-functionEvaluation + 1) / 2;
}

static inline void function8(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
{
    functionEvaluation = 8 * ((abs(x) - 0.5)).square() - 1;
    slope = 2 * (abs(x) - 0.5) * x.sign();
}

static inline void function9(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument)
{
    if (optionalFunctionArgument - round(optionalFunctionArgument) == 0) {
        functionEvaluation = pow(1 - abs(x), (int) optionalFunctionArgument) * 2 - 1;
        slope = -x * pow(1 - abs(x), (int) optionalFunctionArgument - 1) / abs(x);
    } else {
        functionEvaluation = pow(1 - abs(x), optionalFunctionArgument) * 2 - 1;
        slope = -x * pow(1 - abs(x), optionalFunctionArgument - 1) / abs(x);
    }
}

void InverseSpaceTransform::onePeriodicFunction(ArrayXXf& x)
{
//    cout << x<<endl<<endl;
    x = x - round(x);   //catastrophic cancellation possible
//    cout << x << endl << endl;
    closeToPeak = abs(x) < maxCloseToPeakDeviation;

    switch (functionSelection) {
        case 1:
            function1(x, functionEvaluation, slope, optionalFunctionArgument);
            break;
        case 2:
            function2(x, functionEvaluation, slope);
            break;
        case 3:
            function3(x, functionEvaluation, slope);
            break;
        case 4:
            function4(x, functionEvaluation, slope);
            break;
        case 5:
            function5(x, functionEvaluation, slope);
            break;
        case 6:
            function6(x, functionEvaluation, slope);
            break;
        case 7:
            function7(x, functionEvaluation, slope, optionalFunctionArgument);
            break;
        case 8:
            function8(x, functionEvaluation, slope);
            break;
        case 9:
            function9(x, functionEvaluation, slope, optionalFunctionArgument);
            break;
        default:
            stringstream errStream;
            errStream << "Selected function is not available.";
            throw BadInputException(errStream.str());
    }
}

void InverseSpaceTransform::setPointsToTransform(const Matrix3Xf& pointsToTransform)
{
    resultsUpToDate = false;
    this->pointsToTransform = pointsToTransform;

    if (pointsToTransform.cols() != pointsToTransformWeights.cols()) {  // actually not a good choise... 
        pointsToTransformWeights_userPreset = RowVectorXf::Ones(pointsToTransform.cols());
        update_pointsToTransformWeights();
    }
}

void InverseSpaceTransform::setPointsToTransformWeights(const RowVectorXf pointsToTransformWeights)
{
    resultsUpToDate = false;
    pointsToTransformWeights_userPreset = pointsToTransformWeights;
    update_pointsToTransformWeights();
}

void InverseSpaceTransform::update_pointsToTransformWeights()
{
    if (radialWeighting) {
        Array< float, 1, Dynamic > radialWeight = pointsToTransform.colwise().squaredNorm().array().rsqrt();
        pointsToTransformWeights = (pointsToTransformWeights_userPreset.array() * radialWeight).matrix();
    } else {
        pointsToTransformWeights = pointsToTransformWeights_userPreset;
    }
    inverseTransformEvaluationScalingFactor = pointsToTransformWeights.sum();
}

void InverseSpaceTransform::setFunctionSelection(int functionSelection)
{
    if (this->functionSelection != functionSelection) {
        resultsUpToDate = false;
        this->functionSelection = functionSelection;
    }
}
void InverseSpaceTransform::setOptionalFunctionArgument(float optionalFunctionArgument)
{
    if (this->optionalFunctionArgument != optionalFunctionArgument) {
        resultsUpToDate = false;
        this->optionalFunctionArgument = optionalFunctionArgument;
    }
}
void InverseSpaceTransform::setLocalTransformFlag()
{
    if (localTransform == false) {
        resultsUpToDate = false;
        localTransform = true;
    }
}
void InverseSpaceTransform::clearLocalTransformFlag()
{
    if (localTransform == true) {
        resultsUpToDate = false;
        localTransform = false;
    }
}
void InverseSpaceTransform::setRadialWeightingFlag()
{
    if (radialWeighting == false) {
        resultsUpToDate = false;
        radialWeighting = true;
    }
    
    update_pointsToTransformWeights();
}
void InverseSpaceTransform::clearRadialWeightingFlag()
{
    if (radialWeighting == true) {
        resultsUpToDate = false;
        radialWeighting = false;
    }
    
    update_pointsToTransformWeights();
}

void InverseSpaceTransform::setMaxCloseToPeakDeviation(float maxCloseToPeakDeviation)
{
    assert(maxCloseToPeakDeviation < 0.5);
    if (this->maxCloseToPeakDeviation != maxCloseToPeakDeviation) {
        resultsUpToDate = false;
        this->maxCloseToPeakDeviation = maxCloseToPeakDeviation;
    }
}

Eigen::Matrix3Xf& InverseSpaceTransform::getGradient()
{
    if (resultsUpToDate) {
        return gradient;
    } else {
        stringstream errStream;
        errStream << "Gradient not up to date, call performTransform() first.";
        throw BadInputException(errStream.str());
    }
}
Eigen::RowVectorXf& InverseSpaceTransform::getInverseTransformEvaluation()
{
    if (resultsUpToDate) {
        return inverseTransformEvaluation;
    } else {
        stringstream errStream;
        errStream << "InverseTransformEvaluation not up to date, call performTransform() first.";
        throw BadInputException(errStream.str());
    }
}
Eigen::RowVectorXf& InverseSpaceTransform::getCloseToPeaksCount()
{
    if (resultsUpToDate) {
        return closeToPeaksCount;
    } else {
        stringstream errStream;
        errStream << "CloseToPeaksCount not up to date, call performTransform() first.";
        throw BadInputException(errStream.str());
    }
}

vector< vector< uint16_t > >& InverseSpaceTransform::getPeaksCloseToEvaluationPositions_indices()
{
    if (resultsUpToDate) {
        peaksCloseToEvaluationPositions_indices.resize(closeToPeak.cols());

        for (int evaluationPositionIndex = 0; evaluationPositionIndex < closeToPeak.cols(); evaluationPositionIndex++) {
            for (int pointIndex = 0; pointIndex < closeToPeak.rows(); pointIndex++) {
                if (closeToPeak(pointIndex, evaluationPositionIndex)) {
                    peaksCloseToEvaluationPositions_indices[evaluationPositionIndex].push_back(pointIndex);
                }
            }
        }

        return peaksCloseToEvaluationPositions_indices;
    } else {
        stringstream errStream;
        errStream << "closeToPeak not up to date, call performTransform() first.";
        throw BadInputException(errStream.str());
    }
}
