/*
 * samplePointsGenerator.h
 *
 *  Created on: 10.04.2017
 *      Author: Yaro
 */

#ifndef SAMPLEPOINTSGENERATOR_H_
#define SAMPLEPOINTSGENERATOR_H_

#include <string.h>
#include <Eigen/Dense>

class SamplePointsGenerator {
public:
    SamplePointsGenerator();
    SamplePointsGenerator(const std::string& precomputedSamplePointsPath);

    void getDenseGrid(Eigen::Matrix3Xf& samplePoints, float unitPitch, float minRadius, float maxRadius);

    void getTightGrid(Eigen::Matrix3Xf& samplePoints, float unitPitch, float tolerance, const Eigen::VectorXf radii);

private:
    std::string precomputedSamplePointsPath;

    void loadPrecomputedSamplePoints(Eigen::Matrix3Xf& samplePoints, float unitPitch, float tolerance);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* SAMPLEPOINTSGENERATOR_H_ */
