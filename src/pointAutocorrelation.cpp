/*
 * pointAutocorrelation.cpp
 *
 *  Created on: 06.05.2017
 *      Author: Yaro
 */

#include <pointAutocorrelation.h>
#include <ctype.h>

using namespace std;
using namespace Eigen;

void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, const Matrix3Xf& points, float maxNormInAutocorrelation)
{

}

void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, VectorXf& centerPointIndices, VectorXf& shiftedPointIndices,
        const Matrix3Xf& points, float maxNormInAutocorrelation)
{
    uint32_t N = points.cols();
    uint32_t n = N - 1;
    uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

    autocorrelationPoints.resize(3, maxAutocorrPointsCount);
    centerPointIndices.resize(maxAutocorrPointsCount);
    shiftedPointIndices.resize(maxAutocorrPointsCount);

    uint32_t autocorrelationPointsCount = 0;

    for (uint32_t i = 0; i < N; ++i) {
        uint32_t addedNodesInIteration = n - i;

        autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) = points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

        for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k) {
            if (autocorrelationPoints(0, j) > 0) {
                centerPointIndices[j] = i;
                shiftedPointIndices[j] = k;
            } else {
                autocorrelationPoints.col(j) *= -1;
                centerPointIndices[j] = k;
                shiftedPointIndices[j] = i;
            }
        }

        autocorrelationPointsCount += addedNodesInIteration;
    }

    uint32_t goodValuesCount = 0;

    if (maxNormInAutocorrelation == numeric_limits< float >::max()) {
        goodValuesCount = maxAutocorrPointsCount;
    } else {
        float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;

        Array< bool, 1, Eigen::Dynamic > goodValues = autocorrelationPoints.colwise().squaredNorm().array() < maxNormInAutocorrelation_squared;

        for (uint32_t i = 0; i < maxAutocorrPointsCount; ++i) {
            if (goodValues[i]) {
                autocorrelationPoints.col(goodValuesCount) = autocorrelationPoints.col(i);
                centerPointIndices[goodValuesCount] = centerPointIndices[i];
                shiftedPointIndices[goodValuesCount] = shiftedPointIndices[i];
                goodValuesCount++;
            }
        }
    }

    autocorrelationPoints.conservativeResize(3, goodValuesCount);
    centerPointIndices.conservativeResize(goodValuesCount);
    shiftedPointIndices.conservativeResize(goodValuesCount);
}
