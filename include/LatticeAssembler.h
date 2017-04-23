/*
 * LatticeAssembler.h
 *
 *  Created on: 22.04.2017
 *      Author: Yaro
 */

#ifndef LATTICEASSEMBLER_H_
#define LATTICEASSEMBLER_H_

#include "Lattice.h"
#include <vector>
#include <Eigen/Dense>

class LatticeAssembler {
public:
    LatticeAssembler(Eigen::Vector2f& determinantRange, int minNodesOnBasis);

private:
    //input
    Eigen::Vector2f determinantRange;
    int minPointsOnBasis;

    //output
    Eigen::Matrix3Xf pointsToFitInInverseSpace;
    std::vector< Lattice > validLattices;

    //internal
    typedef struct {
        Lattice realSpaceLattice;
        float weight;
        std::vector< uint16_t > pointOnLatticeIndices;
        uint16_t vectorIndices[3];

        uint16_t occupiedLatticePointsCount;
        float meanDefect;
        float meanRelativeDefect;
    } candidateLattice_t;

    std::vector< candidateLattice_t > candidateLattices;
    void computeCandidateLattices(Eigen::Matrix3Xf candidateVectors, Eigen::RowVectorXf candidateVectorWeights,std::vector< std::vector< uint16_t > > pointIndicesOnVector);
    void computeSecondaryLatticeStatistics(candidateLattice_t& candidateLattice);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* LATTICEASSEMBLER_H_ */
