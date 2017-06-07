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
#include <fstream>
#include <eigenSTLContainers.h>
#include "pointAutocorrelation.h"

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
    precomputeIndexingStrategy_balanced();
    precomputeIndexingStrategy_autocorrPrefit();
}

void Indexer::precomputeIndexingStrategy_balanced()
{
    if (experimentSettings.isLatticeParametersKnown()) {
        float unitPitch = 0.05;
        float tolerance = 0.02;

        samplePointsGenerator.getTightGrid(samplePoints_balanced, unitPitch, tolerance, experimentSettings.getDifferentRealLatticeVectorLengths_A());
    } else {
        float unitPitch = 0.05;
        float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
        float maxRadius = experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

        samplePointsGenerator.getDenseGrid(samplePoints_balanced, unitPitch, minRadius, maxRadius);
    }

    float minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.2;
    float maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.2;
    sparsePeakFinder_balanced.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);
}

void Indexer::index_balanced(vector< Lattice >& assembledLattices, const Matrix2Xf& detectorPeaks_m)
{
    if (samplePoints_balanced.size() == 0) {
        precomputeIndexingStrategy_balanced();
    }

    Matrix3Xf samplePoints = samplePoints_balanced;

    Matrix3Xf reciprocalPeaks_A;
    detectorToReciprocalSpaceTransform.computeReciprocalPeaksFromDetectorPeaks(reciprocalPeaks_A, detectorPeaks_m);

    //////// global hill climbing
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
    hillClimbing_accuracyConstants_global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_global);
    hillClimbingOptimizer.performOptimization(reciprocalPeaks_A, samplePoints);

//    ofstream ofs0("workfolder/weights0", ofstream::out);
//    ofs0 << hillClimbingOptimizer.getLastInverseTransformEvaluation().transpose().eval();

//    ofstream ofs("workfolder/samplePoints", ofstream::out);
//    ofs << samplePoints.transpose().eval();

/////// keep only big values 
//    float minFunctionEvaluation = 0.55; //not needed here, since C++ peak finder is much faster than Matlab peak finder.
//    RowVectorXf& samplePointsEvaluation = hillClimbingOptimizer.getLastInverseTransformEvaluation();
//    clearSamplePointsWithLowInverseFunctionEvaluation(samplePoints, samplePointsEvaluation, minFunctionEvaluation);

    /////// find peaks
    uint32_t maxPeaksToTakeCount = 50;
    sparsePeakFinder_balanced.findPeaks_fast(samplePoints, hillClimbingOptimizer.getLastInverseTransformEvaluation());

//    ofstream ofs("workfolder/samplePoints", ofstream::out);
//    ofs << samplePoints.transpose().eval();

    filterSamplePointsForInverseFunctionEvaluation(samplePoints, hillClimbingOptimizer.getLastInverseTransformEvaluation(), maxPeaksToTakeCount);

    /////// peaks hill climbing
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_peaks;

    hillClimbing_accuracyConstants_peaks.functionSelection = 9;
    hillClimbing_accuracyConstants_peaks.optionalFunctionArgument = 8;
    hillClimbing_accuracyConstants_peaks.maxCloseToPeakDeviation = 0.15;

    hillClimbing_accuracyConstants_peaks.initialIterationCount = 0;
    hillClimbing_accuracyConstants_peaks.calmDownIterationCount = 0;
    hillClimbing_accuracyConstants_peaks.calmDownFactor = 0;
    hillClimbing_accuracyConstants_peaks.localFitIterationCount = 10;
    hillClimbing_accuracyConstants_peaks.localCalmDownIterationCount = 20;
    hillClimbing_accuracyConstants_peaks.localCalmDownFactor = 0.85;

    hillClimbing_accuracyConstants_peaks.stepComputationAccuracyConstants.gamma = 0.1;
    hillClimbing_accuracyConstants_peaks.stepComputationAccuracyConstants.maxStep = experimentSettings.getDifferentRealLatticeVectorLengths_A().mean() / 2000;
    hillClimbing_accuracyConstants_peaks.stepComputationAccuracyConstants.minStep = experimentSettings.getDifferentRealLatticeVectorLengths_A().mean() / 20000;
    hillClimbing_accuracyConstants_peaks.stepComputationAccuracyConstants.directionChangeFactor = 2.5;

    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_peaks);
    hillClimbingOptimizer.performOptimization(reciprocalPeaks_A, samplePoints);

//    ofstream ofs("workfolder/samplePoints", ofstream::out);
//    ofs << samplePoints.transpose().eval();

//    ofstream ofs("workfolder/weights1", ofstream::out);
//    ofs << hillClimbingOptimizer.getLastInverseTransformEvaluation().transpose().eval();

    /////// assemble lattices
    latticeAssembler.reset();

    LatticeAssembler::accuracyConstants_t accuracyConstants_LatticeAssembler;
    accuracyConstants_LatticeAssembler.maxCountGlobalPassingWeightFilter = 500;
    accuracyConstants_LatticeAssembler.maxCountLocalPassingWeightFilter = 15;
    accuracyConstants_LatticeAssembler.maxCountPassingRelativeDefectFilter = 50;
    accuracyConstants_LatticeAssembler.minPointsOnLattice = 5;

    if (experimentSettings.isLatticeParametersKnown()) { //TODO: bad design! put tolerance in experiment settings and compute getMinRealLatticeDeterminant_A3 from that!!!
        latticeAssembler.setDeterminantRange(experimentSettings.getRealLatticeDeterminant_A3() * 0.8, experimentSettings.getRealLatticeDeterminant_A3() * 1.2);
    } else {
        latticeAssembler.setDeterminantRange(experimentSettings.getMinRealLatticeDeterminant_A3(), experimentSettings.getMaxRealLatticeDeterminant_A3());
    }

    latticeAssembler.setAccuracyConstants(accuracyConstants_LatticeAssembler);

    vector< LatticeAssembler::assembledLatticeStatistics_t > assembledLatticesStatistics;
    Matrix3Xf& candidateVectors = samplePoints;
    RowVectorXf& candidateVectorWeights = hillClimbingOptimizer.getLastInverseTransformEvaluation(); //TODO: FEHLER!!! u.U. falsche radial weightening und local transform flag!!
    vector< vector< uint16_t > >& pointIndicesOnVector = hillClimbingOptimizer.getPeaksCloseToEvaluationPositions_indices();
    latticeAssembler.assembleLattices(assembledLattices, assembledLatticesStatistics, candidateVectors,
            candidateVectorWeights, pointIndicesOnVector, reciprocalPeaks_A);

//    cout << assembledLatticesStatistics[0].meanDefect << " " << assembledLatticesStatistics[0].meanRelativeDefect << " "
//            << assembledLatticesStatistics[0].occupiedLatticePointsCount << " " << assembledLatticesStatistics.size() << endl << assembledLattices[0].det()
//            << endl;
}

void Indexer::getGoodAutocorrelationPoints(Matrix3Xf& goodAutocorrelationPoints, RowVectorXf& goodAutocorrelationPointWeights, const Matrix3Xf& points,
        uint32_t maxAutocorrelationPointsCount)
{
    Matrix3Xf autocorrelationPoints;

    getPointAutocorrelation(autocorrelationPoints, points, minNormInAutocorrelation_autocorrPrefit, maxNormInAutocorrelation_autocorrPrefit);

    vector< Dbscan::cluster_t > clusters;
    uint16_t minPoints = 2;
    dbscan_autocorrPrefit.computeClusters(clusters, autocorrelationPoints, minPoints, dbscanEpsilon_autocorrPrefit);

    sort(clusters.begin(), clusters.end(), [&](const Dbscan::cluster_t& i, const Dbscan::cluster_t& j) {return i.size() > j.size();});

    Array< uint8_t, 1, Dynamic > autocorrelationPointIsInCluster(autocorrelationPoints.count());
    for (auto& cluster : clusters) {
        for (uint32_t index : cluster) {
            autocorrelationPointIsInCluster[index] = 1;
        }
    }
    uint32_t pointsOutsideOfClustersCount = autocorrelationPointIsInCluster.size() - autocorrelationPointIsInCluster.sum();

    goodAutocorrelationPoints.resize(3, min(pointsOutsideOfClustersCount + (uint32_t) clusters.size(), maxAutocorrelationPointsCount));

    uint32_t goodAutocorrelationPointsCount = 0;

    uint32_t clusterMeansToTakeCount = min((uint32_t) goodAutocorrelationPoints.size(), (uint32_t) clusters.size());
    for (; goodAutocorrelationPointsCount < clusterMeansToTakeCount; ++goodAutocorrelationPointsCount) {
        Vector3f sum(0, 0, 0);
        for (uint32_t index : clusters[goodAutocorrelationPointsCount]) {
            sum += autocorrelationPoints.col(index);
        }
        Vector3f mean = sum / clusters[goodAutocorrelationPointsCount].size();

        goodAutocorrelationPoints.col(goodAutocorrelationPointsCount) = mean;
        goodAutocorrelationPointWeights[goodAutocorrelationPointsCount] = clusters[goodAutocorrelationPointsCount].size();
    }

    uint32_t pointsOutsideOfClustersToTakeCount = goodAutocorrelationPoints.cols() - goodAutocorrelationPointsCount;
    if (pointsOutsideOfClustersToTakeCount == 0) {
        return;
    }

    EigenSTL::vector_Vector3f pointsOutsideOfClusters;
    pointsOutsideOfClusters.reserve(pointsOutsideOfClustersCount);
    for (int i = 0; i < autocorrelationPoints.size(); ++i) {
        if (!autocorrelationPointIsInCluster[i]) {
            pointsOutsideOfClusters.push_back(autocorrelationPoints.col(i));
        }
    }

    nth_element(pointsOutsideOfClusters.begin(), pointsOutsideOfClusters.begin() + pointsOutsideOfClustersToTakeCount, pointsOutsideOfClusters.end(),
            [&](const Vector3f& i, const Vector3f& j) {return i.squaredNorm() < j.squaredNorm();});

    for (uint32_t i = 0; i < pointsOutsideOfClustersToTakeCount; ++i, ++goodAutocorrelationPointsCount) {
        goodAutocorrelationPoints.col(goodAutocorrelationPointsCount) = pointsOutsideOfClusters[i];
        goodAutocorrelationPointWeights[goodAutocorrelationPointsCount] = 1;
    }

    goodAutocorrelationPointWeights.array().cube().sqrt(); //weighten bigger clusters more
}

void Indexer::precomputeIndexingStrategy_autocorrPrefit()
{
    if (experimentSettings.isLatticeParametersKnown()) {
        float unitPitch = 0.05;
        float tolerance = 0.02;

        samplePointsGenerator.getTightGrid(samplePoints_autocorrPrefit, unitPitch, tolerance, experimentSettings.getDifferentRealLatticeVectorLengths_A());
    } else {
        float unitPitch = 0.05;
        float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
        float maxRadius = experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

        samplePointsGenerator.getDenseGrid(samplePoints_autocorrPrefit, unitPitch, minRadius, maxRadius);
    }

    maxNormInAutocorrelation_autocorrPrefit = experimentSettings.getMaxReciprocalLatticeVectorLength_1A() * 5;
    minNormInAutocorrelation_autocorrPrefit = experimentSettings.getMinReciprocalLatticeVectorLength_1A() * 0.7;
    dbscanEpsilon_autocorrPrefit = experimentSettings.getMinReciprocalLatticeVectorLength_1A() * 0.15;
    dbscan_autocorrPrefit.init(dbscanEpsilon_autocorrPrefit, maxNormInAutocorrelation_autocorrPrefit);

    float minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.3;
    float maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.2;
    sparsePeakFinder_autocorrPrefit.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);

    maxCloseToPeakDeviation_autocorrPrefit = 0.15;
    inverseSpaceTransform_autocorrPrefit = InverseSpaceTransform(maxCloseToPeakDeviation_autocorrPrefit);
}

void Indexer::autocorrPrefit(const Matrix3Xf& reciprocalPeaks_A, Matrix3Xf& samplePoints,
        HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_autocorr)
{
    Matrix3Xf autocorrelationPeaks;
    RowVectorXf autocorrelationPointWeights;

    float maxAutocorrelationPointsCount = 20;

    getGoodAutocorrelationPoints(autocorrelationPeaks, autocorrelationPointWeights, reciprocalPeaks_A, maxAutocorrelationPointsCount);

    if (autocorrelationPeaks.cols() < 5) {
        return;
    }

    hillClimbing_accuracyConstants_autocorr.functionSelection = 1;
    hillClimbing_accuracyConstants_autocorr.optionalFunctionArgument = 1;
    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_autocorr);
    hillClimbingOptimizer.performOptimization(autocorrelationPeaks, samplePoints);

    Matrix3Xf peakPositions_11 = samplePoints;
    RowVectorXf peakValues_11 = hillClimbingOptimizer.getLastInverseTransformEvaluation();
    sparsePeakFinder_autocorrPrefit.findPeaks_fast(peakPositions_11, peakValues_11);
}

void Indexer::index_autocorrPrefit(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m)
{
    if (samplePoints_autocorrPrefit.size() == 0) {
        precomputeIndexingStrategy_autocorrPrefit();
    }

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
            [&](uint32_t i, uint32_t j) {return samplePointsEvaluation[i] > samplePointsEvaluation[j];});

    sortIndices.resize(toTakeCount);
    sort(sortIndices.begin(), sortIndices.end(),
            [&](uint32_t i, uint32_t j) {return samplePointsEvaluation[i] > samplePointsEvaluation[j];});

    Matrix3Xf samplePoints_filtered(3, toTakeCount);
    RowVectorXf samplePointsEvaluation_filtered(toTakeCount);
    for (uint32_t i = 0; i < toTakeCount; ++i) {
        samplePoints_filtered.col(i) = samplePoints.col(sortIndices[i]);
        samplePointsEvaluation_filtered[i] = samplePointsEvaluation[sortIndices[i]];
    }

    samplePoints.swap(samplePoints_filtered);
    samplePointsEvaluation.swap(samplePointsEvaluation_filtered);
}
