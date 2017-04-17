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
    InverseSpaceTransform();
    InverseSpaceTransform(float maxCloseToPeakDeviation);

    void performTransform(const Eigen::Matrix3Xf& positionsToEvaluate);

    void setPointsToTransform(const Eigen::Matrix3Xf& pointsToTransform);
    void setPointsToTransformWeights(const Eigen::RowVectorXf pointsToTransformWeights);

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

    std::vector< std::vector< uint16_t > >& getPeaksCloseToEvaluationPositions_indices();

private:
    void onePeriodicFunction(Eigen::ArrayXXf& x);

    Eigen::Matrix3Xf pointsToTransform;
    Eigen::RowVectorXf pointsToTransformWeights;

    

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
    
    std::vector< std::vector< uint16_t > > peaksCloseToEvaluationPositions_indices;

    bool resultsUpToDate;
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ;
};

#endif /* INVERSESPACETRANSFORM_H_ */
