/*
 * pointAutocorrelation.cpp
 *
 *  Created on: 06.05.2017
 *      Author: Yaro
 */

#include <ctype.h>
#include <pointAutocorrelation.h>

using namespace std;
using namespace Eigen;

void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, const Matrix3Xf& points, float minNormInAutocorrelation, float maxNormInAutocorrelation)
{
    uint32_t N = points.cols();
    uint32_t n = N - 1;
    uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

    float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;
    float minNormInAutocorrelation_squared = minNormInAutocorrelation * minNormInAutocorrelation;

    autocorrelationPoints.resize(3, maxAutocorrPointsCount);

    uint32_t autocorrelationPointsCount = 0;

    for (uint32_t i = 0; i < N; ++i)
    {
        uint32_t addedNodesInIteration = n - i;

        autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) = points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

        for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k)
        {
            if (autocorrelationPoints(2, j) < 0)
            {
                autocorrelationPoints.col(j) *= -1;
            }
        }

        uint32_t endInThisIteration = autocorrelationPointsCount + addedNodesInIteration;
        for (uint32_t j = autocorrelationPointsCount, k = 0; j < endInThisIteration; ++j, ++k)
        {
            float squaredNorm = autocorrelationPoints.col(j).squaredNorm();
            if (squaredNorm < maxNormInAutocorrelation_squared && squaredNorm > minNormInAutocorrelation_squared)
            {
                autocorrelationPoints.col(autocorrelationPointsCount) = autocorrelationPoints.col(j);
                autocorrelationPointsCount++;
            }
        }
    }

    autocorrelationPoints.conservativeResize(3, autocorrelationPointsCount);
}

void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, VectorXi& centerPointIndices, VectorXi& shiftedPointIndices, const Matrix3Xf& points,
                             float minNormInAutocorrelation, float maxNormInAutocorrelation)
{
    uint32_t N = points.cols();
    uint32_t n = N - 1;
    uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

    float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;
    float minNormInAutocorrelation_squared = minNormInAutocorrelation * minNormInAutocorrelation;

    autocorrelationPoints.resize(3, maxAutocorrPointsCount);
    centerPointIndices.resize(maxAutocorrPointsCount);
    shiftedPointIndices.resize(maxAutocorrPointsCount);

    uint32_t autocorrelationPointsCount = 0;

    for (uint32_t i = 0; i < N; ++i)
    {
        uint32_t addedNodesInIteration = n - i;

        autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) = points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

        for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k)
        {
            if (autocorrelationPoints(2, j) >= 0)
            {
                centerPointIndices[j] = i;
                shiftedPointIndices[j] = k;
            }
            else
            {
                autocorrelationPoints.col(j) *= -1;
                centerPointIndices[j] = k;
                shiftedPointIndices[j] = i;
            }
        }

        uint32_t endInThisIteration = autocorrelationPointsCount + addedNodesInIteration;
        for (uint32_t j = autocorrelationPointsCount, k = 0; j < endInThisIteration; ++j, ++k)
        {
            float squaredNorm = autocorrelationPoints.col(j).squaredNorm();
            if (squaredNorm < maxNormInAutocorrelation_squared && squaredNorm > minNormInAutocorrelation_squared)
            {
                autocorrelationPoints.col(autocorrelationPointsCount) = autocorrelationPoints.col(j);
                centerPointIndices[autocorrelationPointsCount] = centerPointIndices[j];
                shiftedPointIndices[autocorrelationPointsCount] = shiftedPointIndices[j];
                autocorrelationPointsCount++;
            }
        }
    }

    autocorrelationPoints.conservativeResize(3, autocorrelationPointsCount);
    centerPointIndices.conservativeResize(autocorrelationPointsCount);
    shiftedPointIndices.conservativeResize(autocorrelationPointsCount);
}
