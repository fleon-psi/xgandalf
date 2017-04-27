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
#include <array>
#include <list>

class LatticeAssembler {
public:
    typedef struct {
        uint16_t occupiedLatticePointsCount;
        float meanDefect;
        float meanRelativeDefect;
    } assembledLatticeStatistics_t;

    typedef struct {
        uint32_t maxCountGlobalPassingWeightFilter;
        uint32_t maxCountLocalPassingWeightFilter;
        uint32_t maxCountPassingRelativeDefectFilter;
    } accuracyConstants_t;

    LatticeAssembler(Eigen::Vector2f& determinantRange, int minNodesOnBasis);
    LatticeAssembler(Eigen::Vector2f& determinantRange, int minNodesOnBasis, accuracyConstants_t accuracyConstants);
    void assembleLattices(std::vector< Lattice >& assembledLattices, Eigen::Matrix3Xf& candidateVectors, Eigen::RowVectorXf& candidateVectorWeights,
            std::vector< std::vector< uint16_t > >& pointIndicesOnVector);
    void assembleLattices(std::vector< Lattice >& assembledLattices, std::vector< assembledLatticeStatistics_t >& assembledLatticesStatistics,
            Eigen::Matrix3Xf& candidateVectors, Eigen::RowVectorXf& candidateVectorWeights, std::vector< std::vector< uint16_t > >& pointIndicesOnVector);

private:
    //input
    Eigen::Vector2f determinantRange;
    uint16_t minPointsOnLattice;

    //output
    Eigen::Matrix3Xf pointsToFitInInverseSpace;
    std::vector< Lattice > validLattices;

    //accuracy constants
    accuracyConstants_t accuracyConstants;

    //internal
    typedef struct {
        Lattice realSpaceLattice;
        float weight;
        std::vector< uint16_t > pointOnLatticeIndices;
        std::array< uint16_t, 3 > vectorIndices;

        assembledLatticeStatistics_t assembledLatticeStatistics;

        float det;
    } candidateLattice_t;

    std::vector< candidateLattice_t > candidateLattices;
    void computeCandidateLattices(Eigen::Matrix3Xf candidateVectors, Eigen::RowVectorXf candidateVectorWeights,
            std::vector< std::vector< uint16_t > > pointIndicesOnVector);
    void computeAssembledLatticeStatistics(candidateLattice_t& candidateLattice);
    void selectBestLattices(std::vector< Lattice >& assembledLattices, std::list< candidateLattice_t >& finalCandidateLattices);

    std::vector< uint16_t > sortIndices;  //to avoid frequent reallocation

    uint16_t countUniqueColumns(const Eigen::Matrix3Xf& millerIndices);

    void filterCandidateBasesByWeight(uint32_t maxToTake);
    void filterCandidateBasesByMeanRelativeDefect(uint32_t maxToTake);
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* LATTICEASSEMBLER_H_ */
