/*
 * Dbscan.cpp
 *
 *  Created on: 07.05.2017
 *      Author: Yaro
 */

#include <Dbscan.h>
#include <sstream>

using namespace std;
using namespace Eigen;

Dbscan::Dbscan(float maxEpsilon, float maxPossiblePointNorm) :
        maxEpsilon(maxEpsilon)
{
    neighbourhood.reserve(1000);

    binWidth = maxEpsilon;

    binWidth_reciprocal = 1 / binWidth;

    binsPerDimension = 2 * ceil(maxPossiblePointNorm / binWidth) + 1;
    binCount = binsPerDimension * binsPerDimension * binsPerDimension;
    binCountMinus1 = binCount - 1;
    bin1Position.setConstant(-1.0f * binsPerDimension / 2 * binWidth);
    strides << 1, binsPerDimension, binsPerDimension * binsPerDimension;

    usedBins.reserve(binCount);

    discretizationVolume.resize(binCount);

    const uint32_t typicalMaxPointsPerBin = 100;
    for (auto& bin : discretizationVolume) {
        bin.points.reserve(typicalMaxPointsPerBin);
        bin.pointIndices.reserve(typicalMaxPointsPerBin);
        bin.visited.reserve(typicalMaxPointsPerBin);
        bin.isMemberOfCluster.reserve(typicalMaxPointsPerBin);
        bin.visitedAndMemberOfCluster.reserve(typicalMaxPointsPerBin);
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

    neighbourhood.reserve(27 * typicalMaxPointsPerBin);
}

void Dbscan::computeClusters(vector< cluster_t > clusters, const Matrix3Xf points, uint16_t minPoints, float epsilon)
{
    if (epsilon > maxEpsilon) {
        stringstream errStream;
        errStream << "epsilon must be smaller than maxEpsilon" << endl;
        throw WrongUsageException(errStream.str());
    }

    clusters.reserve(32);   //just for performance

    fillDiscretizationVolume(points);

    for (bin_t& bin : usedBins) {
        for (int i = 0; i < bin.points.size(); ++i) {

        }
    }

}

void Dbscan::expandCluster(uint32_t point, cluster_t& cluster)
{
}

uint32_t Dbscan::regionQuery(uint32_t point)
{
}

void Dbscan::fillDiscretizationVolume(const Eigen::Matrix3Xf& points)
{
    uint32_t pointsCount = points.cols();

    for (auto& bin : discretizationVolume) {
        bin.points.clear();
        bin.pointIndices.clear();
        bin.visited.clear();
        bin.isMemberOfCluster.clear();
        bin.visitedAndMemberOfCluster.clear();
    }

    for (uint32_t i = 0; i < pointsCount; ++i) {
        const Vector3f& point = points.col(i);
        const uint32_t index = getIndex(point);
        bin_t& bin = discretizationVolume[index];

        bin.points.emplace_back(point.x(), point.y(), point.z(), 0);
        bin.pointIndices.push_back(i);
        bin.visited.push_back(false);
        bin.isMemberOfCluster.push_back(false);
        bin.visitedAndMemberOfCluster.push_back(false);

        usedBins.insert(&bin);
    }
}
