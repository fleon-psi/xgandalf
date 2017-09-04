/*
 * LatticeAssembler.h
 *
 *  Created on: 22.04.2017
 *      Author: Yaro
 */

#ifndef LATTICEASSEMBLER_H_
#define LATTICEASSEMBLER_H_

#include "Lattice.h"
#include <Eigen/Dense>
#include <array>
#include <list>
#include <vector>

class LatticeAssembler
{
  public:
    typedef struct
    {
        uint16_t occupiedLatticePointsCount;
        float meanDefect;
        float meanRelativeDefect;
    } assembledLatticeStatistics_t;

    typedef struct
    {
        uint32_t maxCountGlobalPassingWeightFilter;
        uint32_t maxCountLocalPassingWeightFilter;
        uint32_t maxCountPassingRelativeDefectFilter;
        uint16_t minPointsOnLattice;
    } accuracyConstants_t;

    LatticeAssembler();
    LatticeAssembler(const Eigen::Vector2f& determinantRange);
    LatticeAssembler(const Eigen::Vector2f& determinantRange, const Lattice& sampleRealLattice_A, float knownLatticeTolerance);
    LatticeAssembler(const Eigen::Vector2f& determinantRange, const accuracyConstants_t& accuracyConstants);
    LatticeAssembler(const Eigen::Vector2f& determinantRange, const Lattice& sampleRealLattice_A, float knownLatticeTolerance,
                     const accuracyConstants_t& accuracyConstants);

    void setAccuracyConstants(const accuracyConstants_t& accuracyConstants);
    void setDeterminantRange(const Eigen::Vector2f& determinantRange);
    void setDeterminantRange(float min, float max);
    void setKnownLatticeParameters(const Lattice& sampleRealLattice_A, float tolerance);

    void assembleLattices(std::vector<Lattice>& assembledLattices, Eigen::Matrix3Xf& candidateVectors, Eigen::RowVectorXf& candidateVectorWeights,
                          std::vector<std::vector<uint16_t>>& pointIndicesOnVector, Eigen::Matrix3Xf& pointsToFitInReciprocalSpace);
    void assembleLattices(std::vector<Lattice>& assembledLattices, std::vector<assembledLatticeStatistics_t>& assembledLatticesStatistics,
                          Eigen::Matrix3Xf& candidateVectors, Eigen::RowVectorXf& candidateVectorWeights,
                          std::vector<std::vector<uint16_t>>& pointIndicesOnVector, Eigen::Matrix3Xf& pointsToFitInReciprocalSpace);

  private:
    void setStandardValues();
    void reset();
    bool checkLatticeParameters(Lattice& lattice);

    // input
    Eigen::Vector2f determinantRange;
    Eigen::Array<float, 6, 1> knownLatticeParameters;
    Eigen::Array<float, 6, 1> knownLatticeParametersInverse;
    float knownLatticeParametersTolerance;
    bool latticeParametersKnown;

    accuracyConstants_t accuracyConstants;

    // internal
    typedef struct
    {
        Lattice realSpaceLattice;
        float weight;
        std::vector<uint16_t> pointOnLatticeIndices;
        std::array<uint16_t, 3> vectorIndices;

        assembledLatticeStatistics_t assembledLatticeStatistics;

        float det;
    } candidateLattice_t;

    std::vector<candidateLattice_t> candidateLattices;
    void computeCandidateLattices(Eigen::Matrix3Xf& candidateVectors, Eigen::RowVectorXf& candidateVectorWeights,
                                  std::vector<std::vector<uint16_t>>& pointIndicesOnVector);
    void computeAssembledLatticeStatistics(candidateLattice_t& candidateLattice, const Eigen::Matrix3Xf& pointsToFitInReciprocalSpace);
    void selectBestLattices(std::vector<Lattice>& assembledLattices, std::vector<assembledLatticeStatistics_t>& assembledLatticesStatistics,
                            std::list<candidateLattice_t>& finalCandidateLattices);

    std::vector<uint32_t> sortIndices; // to avoid frequent reallocation
    std::vector<Lattice> validLattices;

    uint16_t countUniqueColumns(const Eigen::Matrix3Xf& millerIndices);

    void filterCandidateLatticesByWeight(uint32_t maxToTakeCount);
    void filterCandidateBasesByMeanRelativeDefect(uint32_t maxToTakeCount);

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* LATTICEASSEMBLER_H_ */
