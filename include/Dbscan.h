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

    typedef struct {
        EigenSTL::vector_Vector4f points;  //padded for being fixed size vectorizable
        std::vector< uint32_t > pointIndices;

        //maybe better std::vector< uint8_t >
        std::vector< bool > visited;
        std::vector< bool > isMemberOfCluster;
        std::vector< bool > visitedAndMemberOfCluster;
    } bin_t;

    //maxEpsilon defines performance. Best performance for maxEpsilon == epsilon. Too small maxEpsilon means very big discretizationVolume
    Dbscan(float maxEpsilon, float maxPossiblePointNorm);

    void computeClusters(std::vector< cluster_t > clusters, const Eigen::Matrix3Xf points, uint16_t minPoints, float epsilon);

private:

    void expandCluster(uint32_t point, cluster_t& cluster);
    uint32_t regionQuery(uint32_t point); //returns neighbours count. Visited neighbours that are members of a cluster are not included in neighbourhood, but counted

    // for less reallocation
    std::vector< bin_t* > neighbourhood;

    //stuff for dicretizationVolume
    float maxEpsilon;

    std::vector< bin_t > discretizationVolume;
    std::unordered_set< bin_t* > usedBins;

    int32_t neighbourBinIndexOffsets[27];

    float binWidth, binWidth_reciprocal;
    int binsPerDimension;
    uint32_t binCount, binCountMinus1;
    Eigen::Vector3f bin1Position;
    Eigen::Vector3d strides;

    void fillDiscretizationVolume(const Eigen::Matrix3Xf& points);

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
