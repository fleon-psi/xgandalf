/*
 * Dbscan.h
 *
 *  Created on: 07.05.2017
 *      Author: Yaro
 */

#ifndef DBSCAN_H_
#define DBSCAN_H_

#include <Eigen/Dense>
#include <vector>
#include <ctype.h>
#include <unordered_set>
#include "eigenSTLContainers.h"
#include "WrongUsageException.h"

class Dbscan {
public:
    typedef std::vector< uint32_t > cluster_t;

    //maxEpsilon defines performance. Best performance for maxEpsilon == epsilon. Too small maxEpsilon means very big discretizationVolume
    Dbscan(float maxEpsilon, float maxPossiblePointNorm);

    void computeClusters(std::vector< cluster_t >& clusters, const Eigen::Matrix3Xf& points, uint16_t minPoints, float epsilon);

private:
    typedef struct {
        EigenSTL::vector_Vector4f points;  //padded for being fixed size vectorizable
        std::vector< uint32_t > pointIndices;

        //maybe better std::vector< uint8_t >
        std::vector< bool > visited;
        std::vector< bool > isMemberOfCluster;
        std::vector< bool > visitedAndMemberOfCluster;
    } bin_t;

    class Neighbour {
    public:
        Neighbour(bin_t* bin, uint32_t indexInBin) :
                bin(bin), indexInBin(indexInBin)
        {
        }

        bin_t* bin;
        uint32_t indexInBin;

        inline bool isVisited() const
        {
            return bin->visited[indexInBin];
        }
        inline bool isMemberOfCluster() const
        {
            return bin->isMemberOfCluster[indexInBin];
        }
        inline uint32_t pointIndex() const
        {
            return bin->pointIndices[indexInBin];
        }
        inline Eigen::Vector4f point() const
        {
            return bin->points[indexInBin];
        }
        inline void markVisited()
        {
            bin->visited[indexInBin] = true;
            if (bin->isMemberOfCluster[indexInBin]) {
                bin->visitedAndMemberOfCluster[indexInBin] = true;
            }
        }
        inline void markIsMemberOfCluster()
        {
            bin->isMemberOfCluster[indexInBin] = true;
            if (bin->visited[indexInBin]) {
                bin->visitedAndMemberOfCluster[indexInBin] = true;
            }
        }
    };

    void expandCluster(cluster_t& cluster, std::vector< Neighbour >& neighbourhood, Neighbour& currentPoint, uint16_t minPoints);
    uint32_t regionQuery(std::vector< Neighbour >& nieghbourhood, const Neighbour& currentPoint); //returns neighbours count. Visited neighbours that are members of a cluster are not included in neighbourhood, but counted

    float squaredEpsilon;

    // for less reallocation
    std::vector< Neighbour > mainNeighbourhood;
    std::vector< Neighbour > neighbourNeighbourhood;

    //stuff for dicretizationVolume
    float maxEpsilon;

    std::vector< bin_t > discretizationVolume;
    std::unordered_set< bin_t* > usedBins;
    //    std::vector< uint32_t > discretizationVolumeIndex;

    int32_t neighbourBinIndexOffsets[27];

    float binWidth, binWidth_reciprocal;
    int binsPerDimension;
    uint32_t binCount, binCountMinus1;
    Eigen::Vector3f bin1Position;
    Eigen::Vector3d strides;

    void fillDiscretizationVolume(const Eigen::Matrix3Xf& points);
    void cleanUpDiscretizationVolume();

    //fast, but insecure: if point lies out of scope, it gets relocated somewhere inside the scope. 
    inline uint32_t getIndex(const Eigen::Vector3f& position) const
            {
        Eigen::Vector3d temp = ((position - bin1Position) * binWidth_reciprocal).cast< double >().array().floor().matrix();
        return std::min((uint32_t) (temp[0] + temp.tail(2).dot(strides.tail(2))), binCountMinus1); //min could be avoided for the price of unsafety
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* DBSCAN_H_ */
