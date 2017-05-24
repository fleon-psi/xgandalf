/*
 * Dbscan.cpp
 *
 *  Created on: 07.05.2017
 *      Author: Yaro
 */

#include <Dbscan.h>
#include <sstream>
#include <iostream>

using namespace std;
using namespace Eigen;

Dbscan::Dbscan(float maxEpsilon, float maxPossiblePointNorm) :
        squaredEpsilon(0), maxEpsilon(maxEpsilon)
{
    if (2 * ceil(maxPossiblePointNorm / maxEpsilon) + 2 + 1) {
        cout << "dbscan histogram would take too much memory! Reducing performance in trade of memory!" << endl << endl;
        this->maxEpsilon = maxPossiblePointNorm / 50;
    }

    binWidth = this->maxEpsilon;

    binWidth_reciprocal = 1 / binWidth;

    binsPerDimension = 2 * ceil(maxPossiblePointNorm / binWidth) + 2 + 1; //+2 for one extra border bin, where nothing should be inside.
    binCount = binsPerDimension * binsPerDimension * binsPerDimension;
    binCountMinus1 = binCount - 1;
    bin1Position.setConstant(-1.0f * binsPerDimension / 2 * binWidth);
    strides << 1, binsPerDimension, binsPerDimension * binsPerDimension;

    discretizationVolume.resize(binCount);

    const uint32_t typicalMaxPointsPerBin = 50;
    for (auto& bin : discretizationVolume) {
        bin.reserve(typicalMaxPointsPerBin);
    }

    neighbourBinIndexOffsets[0] = 0;
    int neighboursIndex = 1;
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            for (int z = -1; z <= 1; z++) {
                if (!(x == 0 && y == 0 && z == 0)) {
                    neighbourBinIndexOffsets[neighboursIndex] = x + y * strides.y() + z * strides.z();
                    neighboursIndex++;
                }
            }
        }
    }

    mainNeighbourhood.reserve(27 * typicalMaxPointsPerBin);
    neighbourNeighbourhood.reserve(27 * typicalMaxPointsPerBin);
}

void Dbscan::computeClusters(vector< cluster_t >& clusters, const Matrix3Xf& points, uint16_t minPoints, float epsilon)
{
    if (epsilon > maxEpsilon) {
        stringstream errStream;
        errStream << "epsilon must be smaller than maxEpsilon" << endl;
        throw WrongUsageException(errStream.str());
    }

    squaredEpsilon = epsilon * epsilon;

    clusters.reserve(50);   //just for performance

    fillDiscretizationVolume(points);

    for (bin_t* bin_p : usedBins) {
        bin_t& bin = *bin_p;
        for (uint32_t i = 0; i < bin.size(); ++i) {
            if (bin[i].visited) {
                continue;
            }

            Neighbour currentPoint(bin_p, &bin[i]);

            currentPoint.markVisited();

            uint32_t neighboursCount = regionQuery(mainNeighbourhood, currentPoint);
            if (neighboursCount >= minPoints) {
                clusters.emplace_back();
                clusters.back().reserve(5);
                expandCluster(clusters.back(), mainNeighbourhood, currentPoint, minPoints);
            }
        }
    }

    cleanUpDiscretizationVolume();
}

void Dbscan::expandCluster(cluster_t& cluster, std::vector< Neighbour >& neighbourhood, Neighbour& currentPoint, uint16_t minPoints)
{
    cluster.push_back(currentPoint.pointIndex());
    currentPoint.markIsMemberOfCluster();

    for (auto neighbour_i = neighbourhood.begin(); neighbour_i != neighbourhood.end(); ++neighbour_i) { //cannot be made a for each loop, since neighbourhood size changes
        auto& neighbour = *neighbour_i;
        if (!neighbour.isVisited()) {
            neighbour.markVisited();
            uint32_t neighboursCount = regionQuery(neighbourNeighbourhood, neighbour);
            if (neighboursCount >= minPoints) {
                neighbourhood.insert(neighbourhood.end(), neighbourNeighbourhood.begin(), neighbourNeighbourhood.end());
            }
        }
        if (!neighbour.isMemberOfCluster()) {
            cluster.push_back(neighbour.pointIndex());
            neighbour.markIsMemberOfCluster();
        }
    }
}

uint32_t Dbscan::regionQuery(std::vector< Neighbour >& nieghbourhood, const Neighbour& currentPoint)
{
    uint32_t validNeighboursCount = 0;
    nieghbourhood.clear();

    const Vector4f& currentPointPos = currentPoint.point();
    for (int neighbourBinIndexOffset : neighbourBinIndexOffsets) {
        bin_t& neighbourBin = *(currentPoint.bin + neighbourBinIndexOffset);
        for (uint32_t i = 0; i < neighbourBin.size(); ++i) {
            const Vector4f& neighbourPos = neighbourBin[i].point;
            if ((currentPointPos - neighbourPos).squaredNorm() <= squaredEpsilon) {
                validNeighboursCount++;
                if (!neighbourBin[i].visitedAndMemberOfCluster) {
                    nieghbourhood.emplace_back(&neighbourBin, &neighbourBin[i]);
                }
            }
        }
    }

    return validNeighboursCount;
}

void Dbscan::fillDiscretizationVolume(const Eigen::Matrix3Xf& points)
{
    uint32_t pointsCount = points.cols();

    for (uint32_t i = 0; i < pointsCount; ++i) {
        const Vector3f& point = points.col(i);
        const uint32_t index = getIndex(point);
        bin_t& bin = discretizationVolume[index];

        bin.emplace_back();
        auto& entry = bin.back();
        entry.point = Vector4f(point.x(), point.y(), point.z(), 0);
        entry.pointIndex = i;
        entry.visited = false;
        entry.isMemberOfCluster = false;
        entry.visitedAndMemberOfCluster = false;

        usedBins.insert(&bin);
    }
}

void Dbscan::cleanUpDiscretizationVolume()
{
    for (bin_t* bin_p : usedBins) {
        (*bin_p).clear();
    }
}
