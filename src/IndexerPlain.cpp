/*
 * IndexerPlain.cpp
 *
 *  Created on: 07.06.2017
 *      Author: Yaro
 */

#include <IndexerPlain.h>
#include <algorithm>
#include <vector>

using namespace Eigen;
using namespace std;

IndexerPlain::IndexerPlain(const ExperimentSettings& experimentSettings)
    : IndexerBase(experimentSettings)
{
    precompute();
}

IndexerPlain::IndexerPlain(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath)
    : IndexerBase(experimentSettings, precomputedSamplePointsPath)
{
    precompute();
}

void IndexerPlain::precompute()
{
    maxCloseToPointDeviation = 0.15;

    setSamplingPitch(SamplingPitch::standard);
    setGradientDescentIterationsCount(GradientDescentIterationsCount::standard);

    float minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.2;
    float maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.2;
    sparsePeakFinder.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);

    inverseSpaceTransform = InverseSpaceTransform(maxCloseToPointDeviation);
    inverseSpaceTransform.setFunctionSelection(9);
    inverseSpaceTransform.setOptionalFunctionArgument(8);
    inverseSpaceTransform.setLocalTransformFlag();
    inverseSpaceTransform.clearRadialWeightingFlag();

    accuracyConstants_LatticeAssembler.maxCountGlobalPassingWeightFilter = 500;
    accuracyConstants_LatticeAssembler.maxCountLocalPassingWeightFilter = 15;
    accuracyConstants_LatticeAssembler.maxCountPassingRelativeDefectFilter = 50;
    accuracyConstants_LatticeAssembler.minPointsOnLattice = 5;
    //    latticeAssembler.setDeterminantRange(experimentSettings.getMinRealLatticeDeterminant_A3(), experimentSettings.getMaxRealLatticeDeterminant_A3());
    latticeAssembler.setDeterminantRange(experimentSettings.getRealLatticeDeterminant_A3() * 0.8,
                                         experimentSettings.getRealLatticeDeterminant_A3() * 1.2); // debug
    latticeAssembler.setAccuracyConstants(accuracyConstants_LatticeAssembler);
    latticeAssembler.setKnownLatticeParameters(experimentSettings.getSampleRealLattice_A(), experimentSettings.getTolerance());
}

void IndexerPlain::setSamplingPitch(SamplingPitch samplingPitch)
{
    float unitPitch = 0;
    bool coverSecondaryMillerIndices = false;

    switch (samplingPitch)
    {
        case SamplingPitch::extremelyLoose:
            unitPitch = 0.10;
            break;
        case SamplingPitch::loose:
            unitPitch = 0.075;
            break;
        case SamplingPitch::standard:
            unitPitch = 0.05;
            break;
        case SamplingPitch::dense:
            unitPitch = 0.025;
            break;
        case SamplingPitch::extremelyDense:
            unitPitch = 0.01;
            break;

        case SamplingPitch::standardWithSeondaryMillerIndices:
            unitPitch = 0.05;
            coverSecondaryMillerIndices = true;
            break;
        case SamplingPitch::denseWithSeondaryMillerIndices:
            unitPitch = 0.025;
            coverSecondaryMillerIndices = true;
            break;
        case SamplingPitch::extremelyDenseWithSeondaryMillerIndices:
            unitPitch = 0.01;
            coverSecondaryMillerIndices = true;
            break;
    }

    setSamplingPitch(unitPitch, coverSecondaryMillerIndices);
}

void IndexerPlain::setSamplingPitch(float unitPitch, bool coverSecondaryMillerIndices)
{
    if (experimentSettings.isLatticeParametersKnown())
    {
        float tolerance = min(unitPitch, experimentSettings.getTolerance());

        if (!coverSecondaryMillerIndices)
        {
            samplePointsGenerator.getTightGrid(precomputedSamplePoints, unitPitch, tolerance, experimentSettings.getDifferentRealLatticeVectorLengths_A());
        }
        else
        {
            auto sampleReaalBasis = experimentSettings.getSampleRealLattice_A().getBasis();

            vector<float> radii(6);
            radii[0] = sampleReaalBasis.col(0).norm();
            radii[1] = sampleReaalBasis.col(1).norm();
            radii[2] = sampleReaalBasis.col(2).norm();
            radii[3] = (sampleReaalBasis.col(0) + sampleReaalBasis.col(1)).norm();
            radii[4] = (sampleReaalBasis.col(1) + sampleReaalBasis.col(2)).norm();
            radii[5] = (sampleReaalBasis.col(2) + sampleReaalBasis.col(0)).norm();

            sort(radii.begin(), radii.end());
            auto it = std::unique(radii.begin(), radii.end());
            radii.resize(std::distance(radii.begin(), it));

            ArrayXf radii_array = Eigen::Map<ArrayXf>(radii.data(), radii.size(), 1);

            samplePointsGenerator.getTightGrid(precomputedSamplePoints, unitPitch, tolerance, radii_array);
        }
    }
    else
    {
        if (!coverSecondaryMillerIndices)
        {
            float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
            float maxRadius = experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

            samplePointsGenerator.getDenseGrid(precomputedSamplePoints, unitPitch, minRadius, maxRadius);
        }
        else
        {
            float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
            float maxRadius = 2 * experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

            samplePointsGenerator.getDenseGrid(precomputedSamplePoints, unitPitch, minRadius, maxRadius);
        }
    }
}

void IndexerPlain::index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A)
{
    if (precomputedSamplePoints.size() == 0)
    {
        precompute();
    }

    Matrix3Xf samplePoints = precomputedSamplePoints;

    // global hill climbing
    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_global);
    hillClimbingOptimizer.performOptimization(reciprocalPeaks_1_per_A, samplePoints);
    RowVectorXf globalHillClimbingPointEvaluation = hillClimbingOptimizer.getLastInverseTransformEvaluation();
    Matrix3Xf globalHillClimbingSamplePoints = samplePoints;

    // additional global hill climbing
    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_additionalGlobal);
    hillClimbingOptimizer.performOptimization(reciprocalPeaks_1_per_A, samplePoints);
    RowVectorXf& additionalGlobalHillClimbingPointEvaluation = hillClimbingOptimizer.getLastInverseTransformEvaluation();
    Matrix3Xf& additionalGlobalHillClimbingSamplePoints = samplePoints;

    // find peaks
    uint32_t maxGlobalPeaksToTakeCount = 50;
    sparsePeakFinder.findPeaks_fast(globalHillClimbingSamplePoints, globalHillClimbingPointEvaluation);
    keepSamplePointsWithHighestEvaluation(globalHillClimbingSamplePoints, globalHillClimbingPointEvaluation, maxGlobalPeaksToTakeCount);

    uint32_t maxAdditionalGlobalPeaksToTakeCount = 50;
    sparsePeakFinder.findPeaks_fast(additionalGlobalHillClimbingSamplePoints, additionalGlobalHillClimbingPointEvaluation);
    keepSamplePointsWithHighestEvaluation(additionalGlobalHillClimbingSamplePoints, additionalGlobalHillClimbingPointEvaluation,
                                          maxAdditionalGlobalPeaksToTakeCount);

    Matrix3Xf peakSamplePoints(3, globalHillClimbingSamplePoints.cols() + additionalGlobalHillClimbingSamplePoints.cols());
    peakSamplePoints << globalHillClimbingSamplePoints, additionalGlobalHillClimbingSamplePoints;

    // peaks hill climbing
    hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_peaks);
    hillClimbingOptimizer.performOptimization(reciprocalPeaks_1_per_A, peakSamplePoints);

    // final peaks extra evaluation
    inverseSpaceTransform.setPointsToTransform(reciprocalPeaks_1_per_A);
    inverseSpaceTransform.performTransform(peakSamplePoints);

    // find peaks , TODO: check, whether better performance without peak finding here
    // sparsePeakFinder.findPeaks_fast(peakSamplePoints, inverseSpaceTransform.getInverseTransformEvaluation());

    // assemble lattices
    vector<LatticeAssembler::assembledLatticeStatistics_t> assembledLatticesStatistics;
    Matrix3Xf& candidateVectors = peakSamplePoints;
    RowVectorXf& candidateVectorWeights = inverseSpaceTransform.getInverseTransformEvaluation();
    vector<vector<uint16_t>>& pointIndicesOnVector = inverseSpaceTransform.getPointsCloseToEvaluationPositions_indices();
    Matrix3Xf reciprocalPeaksCopy_1_per_A = reciprocalPeaks_1_per_A;
    latticeAssembler.assembleLattices(assembledLattices, assembledLatticesStatistics, candidateVectors, candidateVectorWeights, pointIndicesOnVector,
                                      reciprocalPeaksCopy_1_per_A);

    //    cout << assembledLatticesStatistics[0].meanDefect << " " << assembledLatticesStatistics[0].meanRelativeDefect << " "
    //            << assembledLatticesStatistics[0].occupiedLatticePointsCount << " " << assembledLatticesStatistics.size() << endl <<
    //            assembledLattices[0].det()
    //            << endl;

    //    ofstream ofs("workfolder/samplePoints", ofstream::out);
    //    ofs << samplePoints.transpose().eval();
}

void IndexerPlain::setGradientDescentIterationsCount(GradientDescentIterationsCount gradientDescentIterationsCount)
{
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t& global = hillClimbing_accuracyConstants_global;
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t& additionalGlobal = hillClimbing_accuracyConstants_additionalGlobal;
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t& peaks = hillClimbing_accuracyConstants_peaks;

    float meanRealLatticeVectorLength = experimentSettings.getDifferentRealLatticeVectorLengths_A().mean();

    switch (gradientDescentIterationsCount)
    {
        case GradientDescentIterationsCount::exremelyFew:
            global.initialIterationCount = 10;
            global.calmDownIterationCount = 3;
            global.calmDownFactor = 0.65;
            global.localFitIterationCount = 4;
            global.localCalmDownIterationCount = 3;
            global.localCalmDownFactor = 0.65;

            global.stepComputationAccuracyConstants.gamma = 0.75;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 15;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 150;
            global.stepComputationAccuracyConstants.directionChangeFactor = 2;
            break;
        case GradientDescentIterationsCount::few:
            global.initialIterationCount = 20;
            global.calmDownIterationCount = 4;
            global.calmDownFactor = 0.73;
            global.localFitIterationCount = 6;
            global.localCalmDownIterationCount = 5;
            global.localCalmDownFactor = 0.73;

            global.stepComputationAccuracyConstants.gamma = 0.7;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 18;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 180;
            global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

            break;
        case GradientDescentIterationsCount::standard:
            global.initialIterationCount = 40;
            global.calmDownIterationCount = 5;
            global.calmDownFactor = 0.8;
            global.localFitIterationCount = 8;
            global.localCalmDownIterationCount = 6;
            global.localCalmDownFactor = 0.8;

            global.stepComputationAccuracyConstants.gamma = 0.65;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 20;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 200;
            global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

            break;
        case GradientDescentIterationsCount::many:
            global.initialIterationCount = 60;
            global.calmDownIterationCount = 8;
            global.calmDownFactor = 0.83;
            global.localFitIterationCount = 10;
            global.localCalmDownIterationCount = 9;
            global.localCalmDownFactor = 0.83;

            global.stepComputationAccuracyConstants.gamma = 0.65;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 20;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 200;
            global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

            break;
        case GradientDescentIterationsCount::manyMany:
            global.initialIterationCount = 90;
            global.calmDownIterationCount = 10;
            global.calmDownFactor = 0.85;
            global.localFitIterationCount = 15;
            global.localCalmDownIterationCount = 15;
            global.localCalmDownFactor = 0.9;

            global.stepComputationAccuracyConstants.gamma = 0.60;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 23;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 230;
            global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

            break;
        case GradientDescentIterationsCount::extremelyMany:
            global.initialIterationCount = 200;
            global.calmDownIterationCount = 15;
            global.calmDownFactor = 0.9;
            global.localFitIterationCount = 15;
            global.localCalmDownIterationCount = 15;
            global.localCalmDownFactor = 0.9;

            global.stepComputationAccuracyConstants.gamma = 0.60;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 23;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 230;
            global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

            break;
        case GradientDescentIterationsCount::custom:
            global.initialIterationCount = 40;
            global.calmDownIterationCount = 5;
            global.calmDownFactor = 0.8;
            global.localFitIterationCount = 8;
            global.localCalmDownIterationCount = 6;
            global.localCalmDownFactor = 0.8;

            global.stepComputationAccuracyConstants.gamma = 0.65;
            global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 20;
            global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 200;
            global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;
            break;
    }

    global.functionSelection = 1;
    global.optionalFunctionArgument = 1;
    global.maxCloseToPointDeviation = maxCloseToPointDeviation;

    additionalGlobal.functionSelection = 9;
    additionalGlobal.optionalFunctionArgument = 8;
    additionalGlobal.maxCloseToPointDeviation = maxCloseToPointDeviation;

    additionalGlobal.initialIterationCount = 0;
    additionalGlobal.calmDownIterationCount = 0;
    additionalGlobal.calmDownFactor = 0;
    additionalGlobal.localFitIterationCount = 4;
    additionalGlobal.localCalmDownIterationCount = 3;
    additionalGlobal.localCalmDownFactor = 0.7;

    additionalGlobal.stepComputationAccuracyConstants.gamma = 0.65;
    additionalGlobal.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 50;
    additionalGlobal.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 20000;
    additionalGlobal.stepComputationAccuracyConstants.directionChangeFactor = 2.5;

    peaks.functionSelection = 9;
    peaks.optionalFunctionArgument = 8;
    peaks.maxCloseToPointDeviation = maxCloseToPointDeviation;

    peaks.initialIterationCount = 0;
    peaks.calmDownIterationCount = 0;
    peaks.calmDownFactor = 0;
    peaks.localFitIterationCount = 10;
    peaks.localCalmDownIterationCount = 30;
    peaks.localCalmDownFactor = 0.9;

    peaks.stepComputationAccuracyConstants.gamma = 0.1;
    peaks.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 2000;
    peaks.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 20000;
    peaks.stepComputationAccuracyConstants.directionChangeFactor = 2.5;
}
