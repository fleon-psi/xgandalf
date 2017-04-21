/*
 * SparsePeakFinder.h
 *
 *  Created on: 21.04.2017
 *      Author: Yaro
 */

#ifndef SPARSEPEAKFINDER_H_
#define SPARSEPEAKFINDER_H_

#include <ctype.h>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>

class SparsePeakFinder {
public:
    SparsePeakFinder(float minDistanceBetweenRealPeaks, float maxPossiblePointNorm);

    //the only warranty this function makes is, that peaks that are separated by more than minDistanceBetweenRealPeaks are found. Additionally some peaks might be found that are not real peaks
    void findPeaks_weak(Eigen::Matrix3Xf& peakPositions, const Eigen::Matrix3Xf& pointPositions, const Eigen::RowVectorXf& pointValues);

private:
    //fast, but insecure: if point lies out of scope, it gets relocated somewhere inside the scope. Maximum bins per direction is 200;
    inline uint32_t getIndex(Eigen::Vector3f position)
    {
        return std::min((uint32_t) (position - bin1Position).array().round().matrix().dot(strides), binCountMinus1); //implementation with floor might be faster (bin1Position would have to move by half a bin with)
    }

    typedef struct {
        int pointIndex;
        float value;
    } bin_t;

    float minDistanceBetweenRealPeaks;
    float maxPossiblePointNorm;

    float binWidth, reciprocalBinWidth;
    int binsPerDimension;
    uint32_t binCountMinus1;
    Eigen::Vector3f bin1Position;
    Eigen::Vector3f strides;

    std::vector< bin_t > discretizationVolume;
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* SPARSEPEAKFINDER_H_ */
