/*
 * InverseSpaceTransform.h
 *
 *  Created on: 11.04.2017
 *      Author: Yaro
 */

#ifndef INVERSESPACETRANSFORM_H_
#define INVERSESPACETRANSFORM_H_

#include <Eigen/Dense>
#include <ctype.h>
#include <vector>
#include "BadInputException.h"

class InverseSpaceTransform {
public:
    typedef struct {
        int functionSelection;
        float optionalFunctionArgument;
        bool localTransform;
        bool radialWeighting;

        float maxCloseToPeakDeviation;
    } accuracyConstants_t;

    InverseSpaceTransform();
    InverseSpaceTransform(float maxCloseToPeakDeviation);

    void performTransform(const Eigen::Matrix3Xf& positionsToEvaluate);

    void setPointsToTransform(const Eigen::Matrix3Xf& pointsToTransform);
    void setPointsToTransformWeights(const Eigen::RowVectorXf pointsToTransformWeights);

    void setMaxCloseToPeakDeviation(float maxCloseToPeakDeviation);
    void setFunctionSelection(int functionSelection);
    void setOptionalFunctionArgument(float optionalFunctionArgument);
    void setLocalTransformFlag();
    void clearLocalTransformFlag();
    void setRadialWeightingFlag();
    void clearRadialWeightingFlag();

    Eigen::Matrix3Xf& getGradient();
    Eigen::RowVectorXf& getInverseTransformEvaluation();
    Eigen::RowVectorXf& getCloseToPointsCount();

    std::vector< std::vector< uint16_t > >& getPointsCloseToEvaluationPositions_indices();

private:
    void onePeriodicFunction(Eigen::ArrayXXf& x);

    void update_pointsToTransformWeights();

    Eigen::Matrix3Xf pointsToTransform;
    Eigen::RowVectorXf pointsToTransformWeights_userPreset;
    Eigen::RowVectorXf pointsToTransformWeights;

    accuracyConstants_t accuracyConstants;

    //output 
    Eigen::Matrix3Xf gradient;
    Eigen::RowVectorXf inverseTransformEvaluation;
    Eigen::RowVectorXf closeToPointsCount;

    // interna
    Eigen::ArrayXXf functionEvaluation;
    Eigen::ArrayXXf slope;
    Eigen::Array< bool, Eigen::Dynamic, Eigen::Dynamic > closeToPoint;

    std::vector< std::vector< uint16_t > > pointsCloseToEvaluationPositions_indices;

    float inverseTransformEvaluationScalingFactor;

    bool resultsUpToDate;
};

#endif /* INVERSESPACETRANSFORM_H_ */
