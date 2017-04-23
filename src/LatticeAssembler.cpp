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

using namespace Eigen;
using namespace std;

LatticeAssembler::LatticeAssembler(Vector2f& determinantRange, int minNodesOnBasis) :
        determinantRange(determinantRange), minPointsOnBasis(minNodesOnBasis)
{

}

static inline uint16_t countUniqueColumns(const Matrix3Xf& millerIndices)
{
    uint16_t count = 0;

    uint16_t sortIndices[millerIndices.cols()]; //VLA
    iota(sortIndices, sortIndices + millerIndices.cols(), 0);
    sort(sortIndices, sortIndices + millerIndices.cols(), [&millerIndices](int i, int j) {return millerIndices(0,i) < millerIndices(0,j);});

    Matrix3Xf millerIndices_sorted(3, millerIndices.cols() + 1);
    for (int i = 0; i < millerIndices.cols(); ++i) {
        millerIndices_sorted.col(i) = millerIndices.col(sortIndices[i]);
    }
    millerIndices_sorted(0, millerIndices.cols()) = 0.5f; //not a valid miller index! This saves a check in a later loop

    for (int i = 0; i < millerIndices.cols(); ++i) {
        bool isUnique = true;
        for (int j = i + 1; millerIndices_sorted(0, i) == millerIndices_sorted(0, j); j++) {
            if (millerIndices_sorted(1, i) == millerIndices_sorted(1, j) & millerIndices_sorted(2, i) == millerIndices_sorted(2, j)) {
                isUnique = false;
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
        if (pointIndicesOnVector[i].size() < minPointsOnBasis) {
            pointIndicesOnVector[i] = pointIndicesOnVector.back();
            pointIndicesOnVector.pop_back();

            candidateVectors.col(i) = candidateVectors.col(candidateVectors.cols() - 1);
            candidateVectors.conservativeResize(NoChange, candidateVectors.cols() - 1);

            candidateVectorWeights[i] = candidateVectorWeights[candidateVectorWeights.size() - 1];
            candidateVectorWeights.conservativeResize(candidateVectorWeights.size() - 1);
        }
    }

    for_each(pointIndicesOnVector.begin(), pointIndicesOnVector.end(), [](const vector< uint16_t >& v) {sort(v.begin(),v.end());}); //needed for easier intersection computation later

    int candidateVectorsCount = candidateVectors.cols();
    for (int i = 0; i < candidateVectorsCount - 2; ++i) {
        for (int j = (i + 1); j < candidateVectorsCount - 1; ++j) {
            for (int k = (j + 1); k < candidateVectorsCount; ++k) {
                Lattice latticeToCheck(candidateVectors.col(i), candidateVectors.col(j), candidateVectors.col(k));

                float det = latticeToCheck.det();
                if (abs(det) < determinantRange[0] | abs(det) > determinantRange[1]) {
                    continue;
                }

            }
        }
    }
}

void LatticeAssembler::computeSecondaryLatticeStatistics(candidateLattice_t& candidateLattice)
{
//    auto reciprocalBasis = candidateLattice.realSpaceLattice.getBasis().transpose().inverse();
//    Matrix3Xf factorsToReachPoints = reciprocalBasis.inverse()*pointsToFitInInverseSpace;
    Matrix3Xf factorsToReachPoints = candidateLattice.realSpaceLattice.getBasis() * pointsToFitInInverseSpace;

    Matrix3Xf millerIndices = factorsToReachPoints.array().round();

    candidateLattice.occupiedLatticePointsCount = countUniqueColumns(millerIndices);

    auto predictedPoints = candidateLattice.realSpaceLattice.getBasis().inverse() * millerIndices;
    candidateLattice.meanDefect = ((predictedPoints - pointsToFitInInverseSpace).colwise().norm()).mean();
    candidateLattice.meanRelativeDefect = ((factorsToReachPoints - millerIndices).colwise().norm()).mean();
}
