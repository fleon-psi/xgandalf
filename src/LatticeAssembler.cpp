/*
 * LatticeAssembler.cpp
 *
 *  Created on: 22.04.2017
 *      Author: Yaro
 */

#include <LatticeAssembler.h>
#include <ctype.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>

using namespace Eigen;
using namespace std;

LatticeAssembler::LatticeAssembler(Vector2f& determinantRange, int minNodesOnBasis) :
        determinantRange(determinantRange), minPointsOnLattice(minNodesOnBasis)
{
    accuracyConstants.maxCountGlobalPassingWeightFilter = 500;
    accuracyConstants.maxCountLocalPassingWeightFilter = 15;
    accuracyConstants.maxCountPassingRelativeDefectFilter = 50;
}

LatticeAssembler::LatticeAssembler(Vector2f& determinantRange, int minNodesOnBasis, accuracyConstants_t accuracyConstants) :
        determinantRange(determinantRange), minPointsOnLattice(minNodesOnBasis), accuracyConstants(accuracyConstants)
{
}

void LatticeAssembler::assembleLattices(vector< Lattice >& assembledLattices, Matrix3Xf& candidateVectors, RowVectorXf& candidateVectorWeights,
        vector< vector< uint16_t > >& pointIndicesOnVector)
{
    vector< assembledLatticeStatistics_t > assembledLatticesStatistics;
    assembleLattices(assembledLattices, assembledLatticesStatistics, candidateVectors, candidateVectorWeights, pointIndicesOnVector);
}

void LatticeAssembler::assembleLattices(vector< Lattice >& assembledLattices, vector< assembledLatticeStatistics_t >& assembledLatticesStatistics,
        Matrix3Xf& candidateVectors, RowVectorXf& candidateVectorWeights, vector< vector< uint16_t > >& pointIndicesOnVector)
{
    computeCandidateLattices(candidateVectors, candidateVectorWeights, pointIndicesOnVector);

    list< candidateLattice_t > finalCandidateLattices;

    filterCandidateBasesByWeight(accuracyConstants.maxCountGlobalPassingWeightFilter);

    for_each(candidateLattices.begin(), candidateLattices.end(),
            [this](candidateLattice_t& candidateLattice) {
                candidateLattice.realSpaceLattice.minimize();
                candidateLattice.det = candidateLattice.realSpaceLattice.det();
                computeAssembledLatticeStatistics(candidateLattice);
            });

    finalCandidateLattices.insert(finalCandidateLattices.end(), candidateLattices.begin(),
            candidateLattices.begin() + accuracyConstants.maxCountLocalPassingWeightFilter); // assume that candidateVectors is sorted descending for weight!

    filterCandidateBasesByMeanRelativeDefect(accuracyConstants.maxCountPassingRelativeDefectFilter);

    finalCandidateLattices.insert(finalCandidateLattices.end(), candidateLattices.begin(),
            candidateLattices.end());

    selectBestLattices(assembledLattices, finalCandidateLattices);
}

void LatticeAssembler::selectBestLattices(vector< Lattice >& assembledLattices, list< candidateLattice_t >& finalCandidateLattices)
{
    auto& candidateLattices = this->candidateLattices;  //needed for lambda function

    sortIndices.resize(candidateLattices.size());
    iota(sortIndices.begin(), sortIndices.end(), 0);
    sort(sortIndices.begin(), sortIndices.end(),
            [&candidateLattices](uint16_t i, uint16_t j) {return candidateLattices[i].pointOnLatticeIndices.size() > candidateLattices[j].pointOnLatticeIndices.size();}); //descending

    float significantDetReductionFactor = 0.75;
    float significantPointCountReductionFactor = 0.85;
    float significantMeanDefectReductionFactor = 0.7;
    float significantMeanRelativeDefectReductionFactor = 0.8;

    vector< uint16_t > combinedPointIndicesOnSelectedLattices;
    vector< uint16_t > tmp_indices;
    tmp_indices.reserve(1000); //just for speed

    auto bestCandidateLattice = finalCandidateLattices.begin();
    auto nextCandidateLattice = finalCandidateLattices.begin();
    while (bestCandidateLattice != prev(finalCandidateLattices.end()) && bestCandidateLattice != finalCandidateLattices.end()) {
        tmp_indices.clear();
        set_union(combinedPointIndicesOnSelectedLattices.begin(), combinedPointIndicesOnSelectedLattices.end(),
                bestCandidateLattice->pointOnLatticeIndices.begin(), bestCandidateLattice->pointOnLatticeIndices.end(),
                back_inserter(tmp_indices));   // first resizing to take all elements and then rezizing to right size may be slightly faster than back-inserter 
        combinedPointIndicesOnSelectedLattices.swap(tmp_indices);

        nextCandidateLattice = next(bestCandidateLattice);
        while (nextCandidateLattice != finalCandidateLattices.end()) {
            tmp_indices.clear();
            set_difference(nextCandidateLattice->pointOnLatticeIndices.begin(), nextCandidateLattice->pointOnLatticeIndices.end(),
                    combinedPointIndicesOnSelectedLattices.begin(), combinedPointIndicesOnSelectedLattices.end(),
                    back_inserter(tmp_indices));
            uint16_t uniquelyReachedNodesCount = tmp_indices.size();
            if (uniquelyReachedNodesCount >= minPointsOnLattice) { //enough new points on lattice => lattice cannot be rejected
                ++nextCandidateLattice;
            } else if (nextCandidateLattice->pointOnLatticeIndices.size()
                    > bestCandidateLattice->pointOnLatticeIndices.size() * significantPointCountReductionFactor) { //subset of previous lattice + zero to few points; Not significantly fewer points than current lattice
                if (
                (nextCandidateLattice->det <= bestCandidateLattice->det * significantDetReductionFactor
                        && nextCandidateLattice->assembledLatticeStatistics.meanDefect
                                * min(significantMeanDefectReductionFactor, 1.5f / (bestCandidateLattice->det / nextCandidateLattice->det))
                                < bestCandidateLattice->assembledLatticeStatistics.meanDefect
                        && nextCandidateLattice->assembledLatticeStatistics.meanRelativeDefect * significantMeanRelativeDefectReductionFactor
                                < bestCandidateLattice->assembledLatticeStatistics.meanRelativeDefect)
                        ||
                        (nextCandidateLattice->det * significantDetReductionFactor <= bestCandidateLattice->det
                                && nextCandidateLattice->assembledLatticeStatistics.meanDefect < bestCandidateLattice->assembledLatticeStatistics.meanDefect)
                        ||
                        (nextCandidateLattice->assembledLatticeStatistics.meanDefect
                                < bestCandidateLattice->assembledLatticeStatistics.meanDefect
                                        * min(significantMeanDefectReductionFactor, 1.5f / (nextCandidateLattice->det / bestCandidateLattice->det))
                                && nextCandidateLattice->assembledLatticeStatistics.meanRelativeDefect
                                        < bestCandidateLattice->assembledLatticeStatistics.meanRelativeDefect * significantMeanRelativeDefectReductionFactor)
                        ) {
                    auto temp = next(nextCandidateLattice);
                    finalCandidateLattices.splice(bestCandidateLattice, finalCandidateLattices, nextCandidateLattice);
                    nextCandidateLattice = temp;
                } else {
                    nextCandidateLattice = finalCandidateLattices.erase(nextCandidateLattice);
                }
            } else {
                nextCandidateLattice = finalCandidateLattices.erase(nextCandidateLattice);
            }
        }
        ++bestCandidateLattice;
    }

    assembledLattices.clear();
    assembledLattices.reserve(finalCandidateLattices.size());
    for (auto& selectedCandidateLattice : finalCandidateLattices) {
        assembledLattices.push_back(selectedCandidateLattice.realSpaceLattice);
    }

}

uint16_t LatticeAssembler::countUniqueColumns(const Matrix3Xf& millerIndices)
{
    uint16_t count = 0;

    sortIndices.resize(millerIndices.cols());
    iota(sortIndices.begin(), sortIndices.end(), 0);
    sort(sortIndices.begin(), sortIndices.end(), [&millerIndices](uint16_t i, uint16_t j) {return millerIndices(0,i) < millerIndices(0,j);});

    Matrix3Xf millerIndices_sorted(3, millerIndices.cols() + 1);
    for (int i = 0; i < millerIndices.cols(); ++i) {
        millerIndices_sorted.col(i) = millerIndices.col(sortIndices[i]);
    }
    millerIndices_sorted(0, millerIndices.cols()) = 0.5f; //not a valid miller index! This saves an out of bounds check in a later loop

    for (int i = 0; i < millerIndices.cols(); ++i) {
        bool isUnique = true;
        for (int j = i + 1; millerIndices_sorted(0, i) == millerIndices_sorted(0, j); j++) {
            if ((millerIndices_sorted(1, i) == millerIndices_sorted(1, j)) & (millerIndices_sorted(2, i) == millerIndices_sorted(2, j))) {
                isUnique = false;
                break;
            }
        }
        if (isUnique) {
            count++;
        }
    }

    return count;
}

void LatticeAssembler::computeCandidateLattices(Matrix3Xf candidateVectors, RowVectorXf candidateVectorWeights,
        vector< vector< uint16_t > > pointIndicesOnVector)
{
//hand-crafted remove-if for three variables. Delete
    for (int i = candidateVectors.cols() - 1; i >= 0; i--) {
        if (pointIndicesOnVector[i].size() < minPointsOnLattice) {
            pointIndicesOnVector[i] = pointIndicesOnVector.back();
            pointIndicesOnVector.pop_back();

            candidateVectors.col(i) = candidateVectors.col(candidateVectors.cols() - 1);
            candidateVectors.conservativeResize(NoChange, candidateVectors.cols() - 1); //possibly slow, because of call to realloc. can be avoided if needed

            candidateVectorWeights[i] = candidateVectorWeights[candidateVectorWeights.size() - 1];
            candidateVectorWeights.conservativeResize(candidateVectorWeights.size() - 1);
        }
    }

    for_each(pointIndicesOnVector.begin(), pointIndicesOnVector.end(), [](vector< uint16_t >& v) {sort(v.begin(),v.end());}); //needed for easier intersection computation later

    int candidateVectorsCount = candidateVectors.cols();
    for (uint16_t i = 0; i < candidateVectorsCount - 2; ++i) {
        for (uint16_t j = (i + 1); j < candidateVectorsCount - 1; ++j) {
            for (uint16_t k = (j + 1); k < candidateVectorsCount; ++k) {
                Lattice latticeToCheck(candidateVectors.col(i), candidateVectors.col(j), candidateVectors.col(k));

                float absDet = abs(latticeToCheck.det());
                if ((absDet < determinantRange[0]) | (absDet > determinantRange[1])) {
                    continue;
                }

                vector< uint16_t > pointIndicesOnTwoVectors(max(pointIndicesOnVector[i].size(), pointIndicesOnVector[j].size()));
                auto it = set_intersection(pointIndicesOnVector[i].begin(), pointIndicesOnVector[i].end(),
                        pointIndicesOnVector[j].begin(), pointIndicesOnVector[j].end(),
                        pointIndicesOnTwoVectors.begin());
                uint16_t pointsOnBothVectorsCount = it - pointIndicesOnTwoVectors.begin();
                if (pointsOnBothVectorsCount < minPointsOnLattice) {
                    continue;
                }
                pointIndicesOnTwoVectors.resize(pointsOnBothVectorsCount);

                vector< uint16_t > pointIndicesOnLatticeToCheck(max(pointIndicesOnTwoVectors.size(), pointIndicesOnVector[k].size()));
                it = set_intersection(pointIndicesOnTwoVectors.begin(), pointIndicesOnTwoVectors.end(),
                        pointIndicesOnVector[k].begin(), pointIndicesOnVector[k].end(),
                        pointIndicesOnLatticeToCheck.begin());
                uint16_t pointsOnLatticeToCheckCount = it - pointIndicesOnLatticeToCheck.begin();
                if (pointsOnLatticeToCheckCount < minPointsOnLattice) {
                    continue;
                }
                pointIndicesOnLatticeToCheck.resize(pointsOnLatticeToCheckCount);

                candidateLattices.resize(candidateLattices.size() + 1);
                auto& newCandidateBasis = candidateLattices.back();
                newCandidateBasis.realSpaceLattice = latticeToCheck;
                newCandidateBasis.weight = candidateVectorWeights[i] + candidateVectorWeights[j] + candidateVectorWeights[k];
                newCandidateBasis.pointOnLatticeIndices.swap(pointIndicesOnLatticeToCheck);
                newCandidateBasis.vectorIndices = {i,j,k};
            }
        }
    }
}

void LatticeAssembler::computeAssembledLatticeStatistics(candidateLattice_t& candidateLattice)
{
//    auto reciprocalBasis = candidateLattice.realSpaceLattice.getBasis().transpose().inverse();
//    Matrix3Xf factorsToReachPoints = reciprocalBasis.inverse()*pointsToFitInInverseSpace;
    Matrix3Xf factorsToReachPoints = candidateLattice.realSpaceLattice.getBasis() * pointsToFitInInverseSpace;

    Matrix3Xf millerIndices = factorsToReachPoints.array().round();

    candidateLattice.assembledLatticeStatistics.occupiedLatticePointsCount = countUniqueColumns(millerIndices);

    auto predictedPoints = candidateLattice.realSpaceLattice.getBasis().inverse() * millerIndices;
    candidateLattice.assembledLatticeStatistics.meanDefect = ((predictedPoints - pointsToFitInInverseSpace).colwise().norm()).mean();
    candidateLattice.assembledLatticeStatistics.meanRelativeDefect = ((factorsToReachPoints - millerIndices).colwise().norm()).mean();
}

void LatticeAssembler::filterCandidateBasesByWeight(uint32_t maxToTake)
{
    auto& candidateLattices = this->candidateLattices;  //needed for lambda function

    sortIndices.resize(candidateLattices.size());
    iota(sortIndices.begin(), sortIndices.end(), 0);
    nth_element(sortIndices.begin(), sortIndices.begin() + maxToTake, sortIndices.end(),
            [&candidateLattices](uint16_t i, uint16_t j) {return candidateLattices[i].weight > candidateLattices[j].weight;});  //descending
    candidateLattices.resize(min(maxToTake, (uint32_t) candidateLattices.size()));
}

void LatticeAssembler::filterCandidateBasesByMeanRelativeDefect(uint32_t maxToTake)
{
    auto& candidateLattices = this->candidateLattices;  //needed for lambda function

    sortIndices.resize(candidateLattices.size());
    iota(sortIndices.begin(), sortIndices.end(), 0);
    nth_element(sortIndices.begin(), sortIndices.begin() + maxToTake, sortIndices.end(),
            [&candidateLattices](uint16_t i, uint16_t j)
            {   return candidateLattices[i].assembledLatticeStatistics.meanRelativeDefect < candidateLattices[j].assembledLatticeStatistics.meanRelativeDefect;});
    candidateLattices.resize(min(maxToTake, (uint32_t) candidateLattices.size()));
}
