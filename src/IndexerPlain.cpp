/*
 * IndexerPlain.cpp
 *
 *  Created on: 07.06.2017
 *      Author: Yaro
 */

#include <IndexerPlain.h>

using namespace Eigen;
using namespace std;

IndexerPlain::IndexerPlain(const ExperimentSettings& experimentSettings) :
        IndexerBase(experimentSettings)
{
    precompute();
}

IndexerPlain::IndexerPlain(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath) :
        IndexerBase(experimentSettings, precomputedSamplePointsPath)
{
    precompute();
}

//IndexerPlain::~IndexerPlain()
//{
//}

void IndexerPlain::precompute()
{
    if (experimentSettings.isLatticeParametersKnown()) {
        float unitPitch = 0.05;
        float tolerance = 0.02;

        samplePointsGenerator.getTightGrid(samplePoints, unitPitch, tolerance, experimentSettings.getDifferentRealLatticeVectorLengths_A());
    } else {
        float unitPitch = 0.05;
        float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
        float maxRadius = experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

        samplePointsGenerator.getDenseGrid(samplePoints, unitPitch, minRadius, maxRadius);
    }

    float minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.2;
    float maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.2;
    sparsePeakFinder.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);
}

void IndexerPlain::index(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m)
{
    if (samplePoints.size() == 0) {
        precompute();
    }

    Matrix3Xf samplePoints = samplePoints;

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
    sparsePeakFinder.findPeaks_fast(samplePoints, hillClimbingOptimizer.getLastInverseTransformEvaluation());

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
