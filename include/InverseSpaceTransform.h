/*
 * InverseSpaceTransform.h
 *
 *  Created on: 11.04.2017
 *      Author: Yaro
 */

#ifndef INVERSESPACETRANSFORM_H_
#define INVERSESPACETRANSFORM_H_

#include <assert.h>
#include <Eigen/Dense>

class InverseSpaceTransform {
public:
    InverseSpaceTransform();
    InverseSpaceTransform(float maxCloseToPeakDeviation);

    void performTransform(const Eigen::Matrix3Xf& pointsToTransform, const Eigen::Matrix3Xf& positionsToEvaluate);
    void performTransform(const Eigen::Matrix3Xf& pointsToTransform, const Eigen::Matrix3Xf& positionsToEvaluate,
            Eigen::RowVectorXf & pointsToTransformWeights);

    void setFunctionSelection(int functionSelection);
    void setOptionalFunctionArgument(float optionalFunctionArgument);
    void setLocalTransformFlag();
    void clearLocalTransformFlag();
    void setRadialWeightingFlag();
    void clearRadialWeightingFlag();
    void setMaxCloseToPeakDeviation(float maxCloseToPeakDeviation);

    const Eigen::Matrix3Xf getGradient();
    const Eigen::RowVectorXf getInverseTransformEvaluation();
    const Eigen::RowVectorXf getCloseToPeaksCount();
    public:
    void onePeriodicFunction(Eigen::ArrayXXf& x);

    bool resultsUpToDate;

    int functionSelection;
    float optionalFunctionArgument;
    bool localTransform;
    bool radialWeighting;

    float maxCloseToPeakDeviation;

    //output 
    Eigen::Matrix3Xf gradient;
    Eigen::RowVectorXf inverseTransformEvaluation;
    Eigen::RowVectorXf closeToPeaksCount;

    // interna
    Eigen::ArrayXXf functionEvaluation;
    Eigen::ArrayXXf slope;
    Eigen::Array< bool, Eigen::Dynamic, Eigen::Dynamic > closeToPeak;
};

#endif /* INVERSESPACETRANSFORM_H_ */
