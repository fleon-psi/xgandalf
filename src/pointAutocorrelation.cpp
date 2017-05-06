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
    uint32_t N = points.cols();
        uint32_t n = N - 1;
        uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

        float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;

        autocorrelationPoints.resize(3, maxAutocorrPointsCount);

        uint32_t autocorrelationPointsCount = 0;

        for (uint32_t i = 0; i < N; ++i) {
            uint32_t addedNodesInIteration = n - i;

            autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) = points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

            for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k) {
                if (autocorrelationPoints(0, j) < 0) {
                    autocorrelationPoints.col(j) *= -1;
                }
            }

            Array< bool, 1, Eigen::Dynamic > goodValues =
                    autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration).colwise().squaredNorm().array()
                            < maxNormInAutocorrelation_squared;

            uint32_t endInThisIteration = autocorrelationPointsCount + addedNodesInIteration;
            for (uint32_t j = autocorrelationPointsCount, k = 0; j < endInThisIteration; ++j, ++k) {
                if (goodValues[k]) {
                    autocorrelationPoints.col(autocorrelationPointsCount) = autocorrelationPoints.col(j);
                    autocorrelationPointsCount++;
                }
            }
        }

        autocorrelationPoints.conservativeResize(3, autocorrelationPointsCount);
}

void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, VectorXf& centerPointIndices, VectorXf& shiftedPointIndices,
        const Matrix3Xf& points, float maxNormInAutocorrelation)
{
    uint32_t N = points.cols();
    uint32_t n = N - 1;
    uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

    float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;

    autocorrelationPoints.resize(3, maxAutocorrPointsCount);
    centerPointIndices.resize(maxAutocorrPointsCount);
    shiftedPointIndices.resize(maxAutocorrPointsCount);

    uint32_t autocorrelationPointsCount = 0;

    for (uint32_t i = 0; i < N; ++i) {
        uint32_t addedNodesInIteration = n - i;

        autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) = points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

        for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k) {
            if (autocorrelationPoints(0, j) >= 0) {
                centerPointIndices[j] = i;
                shiftedPointIndices[j] = k;
            } else {
                autocorrelationPoints.col(j) *= -1;
                centerPointIndices[j] = k;
                shiftedPointIndices[j] = i;
            }
        }

        Array< bool, 1, Eigen::Dynamic > goodValues =
                autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration).colwise().squaredNorm().array()
                        < maxNormInAutocorrelation_squared;

        uint32_t endInThisIteration = autocorrelationPointsCount + addedNodesInIteration;
        for (uint32_t j = autocorrelationPointsCount, k = 0; j < endInThisIteration; ++j, ++k) {
            if (goodValues[k]) {
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
