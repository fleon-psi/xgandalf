/*
 * SparsePeakFinder.cpp
 *
 *  Created on: 21.04.2017
 *      Author: Yaro
 */

#include <SparsePeakFinder.h>
#include <BadInputException.h>
#include <sstream>
#include <math.h>

using namespace std;
using namespace Eigen;

SparsePeakFinder::SparsePeakFinder(float minDistanceBetweenRealPeaks, float maxPossiblePointNorm)
{
    this->minDistanceBetweenRealPeaks = minDistanceBetweenRealPeaks;
    this->maxPossiblePointNorm = maxPossiblePointNorm;

    binWidth = sqrt(pow(minDistanceBetweenRealPeaks, 2) / 3);   //minDistanceBetweenRealPeaks is diagonal of the cube
    reciprocalBinWidth = 1 / binWidth;

    binsPerDimension = 2 * ceil(maxPossiblePointNorm / binWidth) + 1;
    binCountMinus1 = pow(binsPerDimension, 3) - 1;
    bin1Position.setConstant(-1 * (binsPerDimension - 1) / 2 * binWidth);
    strides << 1, binsPerDimension, pow(binsPerDimension, 2);

    discretizationVolume.resize(binCountMinus1 + 1);

    if (binsPerDimension > 200) {   //pow(2^23,1/3) == 200
        stringstream errStream;
        errStream << "Created discretizationVolume is too big. Function getIndex() will not work (because of the late casting of the float value).";
        throw BadInputException(errStream.str());
    }

}

void SparsePeakFinder::findPeaks_weak(Matrix3Xf& peakPositions, const Matrix3Xf& pointPositions, const RowVectorXf& pointValues)
{
    fill(discretizationVolume.begin(), discretizationVolume.end(), bin_t( { -1, numeric_limits< float >::lowest() }));  //possibly with 0 faster

    for (int pointIndex = 0; pointIndex < pointPositions.cols(); pointIndex++) {
        uint32_t volumeIndex = getIndex(pointPositions.col(pointIndex));
        if (pointValues[pointIndex] > discretizationVolume[volumeIndex].value) {
            discretizationVolume[volumeIndex].value = pointValues[pointIndex];
            discretizationVolume[volumeIndex].pointIndex = pointIndex;
        }
    }
    
    //check for every filled bin whether its point is a peak. Watch out for bins on the border! Possible treating: create one more bin on the border and do not check it. Problem: cannot just go linearly through array.
}
