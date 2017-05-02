/*
 * Indexer.cpp
 *
 *  Created on: 30.04.2017
 *      Author: Yaro
 */

#include <Indexer.h>
#include <ctype.h>
#include <algorithm>
#include <numeric> 

using namespace Eigen;
using namespace std;

Indexer::Indexer(const ExperimentSettings& experimentSettings) :
        experimentSettings(experimentSettings), detectorToReciprocalSpaceTransform(experimentSettings)
{
    construct();
}

Indexer::Indexer(const ExperimentSettings& experimentSettings, const string& precomputedSamplePointsPath) :
        experimentSettings(experimentSettings), samplePointsGenerator(precomputedSamplePointsPath),
                detectorToReciprocalSpaceTransform(experimentSettings)
{
    construct();
}

void Indexer::construct()
{
    precomputeIndexingStrategy_standard();
}

void Indexer::precomputeSamplePoints_standard()
{
    if (experimentSettings.isLatticeParametersKnown()) {
        float unitPitch = 0.05;
        float tolerance = 0.02;

        samplePointsGenerator.getTightGrid(samplePoints_standard, unitPitch, tolerance, experimentSettings.getDifferentRealLatticeVectorLengths_A());
    } else {
        float unitPitch = 0.05;
        float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
        float maxRadius = experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

        samplePointsGenerator.getDenseGrid(samplePoints_standard, unitPitch, minRadius, maxRadius);
    }
}

void Indexer::precomputeIndexingStrategy_standard()
{
    precomputeSamplePoints_standard();

    float minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.2;
    float maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.2;
    sparsePeakFinder_standard.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);
}

void Indexer::index_standard(vector< Lattice >& assembledLattices, const Matrix2Xf& detectorPeaks_m)
{
    if (samplePoints_standard.size() == 0) {
        precomputeSamplePoints_standard();
    }

    Matrix3Xf samplePoints = samplePoints_standard;

    Matrix3Xf reciprocalPeaks_A;
    detectorToReciprocalSpaceTransform.computeReciprocalPeaksFromDetectorPeaks(reciprocalPeaks_A, detectorPeaks_m);

    //////// global hil climbing
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_global;

    hillClimbing_accuracyConstants_global.functionSelection = 1;
    hillClimbing_accuracyConstants_global.optionalFunctionArgument = 1;
    hillClimbing_accuracyConstants_global.maxCloseToPeakDeviation = 0.15;

    hillClimbing_accuracyConstants_global.initialIterationCount = 40;
    hillClimbing_accuracyConstants_global.calmDownIterationCount = 5;
    hillClimbing_accuracyConstants_global.calmDownFactor = 0.8;
    hillClimbing_accuracyConstants_global.localFitIterationCount = 8;
    hillClimbing_accuracyConstants_global.localCalmDownIterationCount = 6;
    hillClimbing_accuracyConstants_global.localCalmDownFactor = 0.8;

    hillClimbing_accuracyConstants_global.stepComputationAccuracyConstants.gamma = 0.65;
    hillClimbing_accuracyConstants_global.stepComputationAccuracyConstants.maxStep = experimentSettings.getDifferentRealLatticeVectorLengths_A().mean() / 20;
    hillClimbing_accuracyConstants_global.stepComputationAccuracyConstants.minStep = experimentSettings.getDifferentRealLatticeVectorLengths_A().mean() / 200;
    hillClimbing_accuracyConstants_global.stepComputationAccuracyConstants.directionChangeFactor = 2.5;

    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_global);
    hillClimbingOptimizer.performOptimization(reciprocalPeaks_A, samplePoints);

    /////// keep only big values 
    float minFunctionEvaluation = 0.65; // TODO: can be much lower, since C++ peak finder is much faster than Matlab peak finder. especially needed for multiple lattices
    RowVectorXf& samplePointsEvaluation = hillClimbingOptimizer.getLastInverseTransformEvaluation();
    clearSamplePointsWithLowInverseFunctionEvaluation(samplePoints, samplePointsEvaluation, minFunctionEvaluation);

    /////// find peaks
    uint32_t maxPeaksToTakeCount = 50;
    sparsePeakFinder_standard.findPeaks_fast(samplePoints, hillClimbingOptimizer.getLastInverseTransformEvaluation());

    filterSamplePointsForInverseFunctionEvaluation(samplePoints, hillClimbingOptimizer.getLastInverseTransformEvaluation(), maxPeaksToTakeCount);
    
    /////// fit peaks
}

void Indexer::clearSamplePointsWithLowInverseFunctionEvaluation(Matrix3Xf& samplePoints, RowVectorXf& samplePointsEvaluation, float minFunctionEvaluation)
{
    uint32_t bigEvaluationSamplePointsCount = 0;
    for (int i = 0; i < samplePointsEvaluation.size(); ++i) {
        if (samplePointsEvaluation[i] >= minFunctionEvaluation) {
            samplePoints.col(bigEvaluationSamplePointsCount) = samplePoints.col(i);
            samplePointsEvaluation[bigEvaluationSamplePointsCount] = samplePointsEvaluation[i];
            bigEvaluationSamplePointsCount++;
        }
    }
    samplePoints.conservativeResize(3, bigEvaluationSamplePointsCount);
    samplePointsEvaluation.conservativeResize(bigEvaluationSamplePointsCount);
}

void Indexer::filterSamplePointsForInverseFunctionEvaluation(Eigen::Matrix3Xf& samplePoints, RowVectorXf& samplePointsEvaluation, uint32_t maxToTakeCount)
{
    uint32_t toTakeCount = min(maxToTakeCount, (uint32_t) samplePointsEvaluation.size());

    sortIndices.resize(samplePointsEvaluation.size());
    iota(sortIndices.begin(), sortIndices.end(), 0);
    nth_element(sortIndices.begin(), sortIndices.begin() + toTakeCount - 1, sortIndices.end(),
            [&](uint32_t i, uint32_t j) {return samplePointsEvaluation[i]> samplePointsEvaluation[j];});

    sortIndices.resize(toTakeCount);
    sort(sortIndices.begin(), sortIndices.end(),
            [&](uint32_t i, uint32_t j) {return samplePointsEvaluation[i]> samplePointsEvaluation[j];});

    Matrix3Xf samplePoints_filtered(3, toTakeCount);
    RowVectorXf samplePointsEvaluation_filtered(toTakeCount);
    for (uint32_t i = 0; i < toTakeCount; ++i) {
        samplePoints_filtered.col(i) = samplePoints.col(i);
        samplePointsEvaluation_filtered[i] = samplePointsEvaluation[sortIndices[i]];
    }

    samplePoints.swap(samplePoints_filtered);
    samplePointsEvaluation.swap(samplePointsEvaluation_filtered);
}
