/*
 * tests.cpp
 *
 *  Created on: 04.05.2017
 *      Author: Yaro
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include "Dbscan.h"
#include "HillClimbingOptimizer.h"
#include "IndexerAutocorrPrefit.h"
#include "IndexerPlain.h"
#include "InverseSpaceTransform.h"
#include "Lattice.h"
#include "LatticeAssembler.h"
#include "SamplePointsGenerator.h"
#include "SparsePeakFinder.h"
#include "eigenDiskImport.h"
#include "pointAutocorrelation.h"
#include "samplePointsFiltering.h"
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "adaptions/crystfel/IndexerPlain.h"

using namespace std;
using namespace Eigen;

static ExperimentSettings getExperimentSettingLys();

void testCrystfelAdaption()
{
    float coffset_m = 0.567855;
    float clen_mm = -439.9992;
    float beamEenergy_eV = 8.0010e+03;
    float divergenceAngle_deg = 0.05 * M_PI / 180;
    float nonMonochromaticity = 0.005;
    float pixelLength_m = 110e-6;
    float detectorRadius_pixel = 750;
    float tolerance = 0.02;

    float detectorDistance_m = clen_mm * 1e-3 + coffset_m;
    float detectorRadius_m = detectorRadius_pixel * pixelLength_m;

    Lattice_t sampleReciprocalLattice_1A = {+0.00945252, -0.00433391, +0.00644485, -0.00298714, -0.01177522,
                                            -0.00374347, +0.01601091, +0.00156280, -0.02065220};

    ExperimentSettings* experimentSettings = ExperimentSettings_new(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg,
                                                                    nonMonochromaticity, sampleReciprocalLattice_1A, tolerance);
    // ExperimentSettings tmp = getExperimentSettingLys();
    // ExperimentSettings* experimentSettings = &tmp;

    char precomputedSamplePointsPath[] = "C:\\DesyFiles\\workspaces\\VisualStudio_workspace\\indexer\\precomputedSamplePoints";
    samplingPitch_t samplingPitch = SAMPLING_PITCH_standard;
    gradientDescentIterationsCount_t gradientDescentIterationsCount = GRADIENT_DESCENT_ITERATION_COUNT_standard;
    IndexerPlain* indexer = IndexerPlain_new(experimentSettings, precomputedSamplePointsPath);
    IndexerPlain_setSamplingPitch(indexer, samplingPitch);
    IndexerPlain_setGradientDescentIterationsCount(indexer, gradientDescentIterationsCount);

    float coordinates_x[20] = {0.027768, 0.02768,  0.0089125, -0.007797, -0.041511, -0.037719,  0.015493, 0.033985, 0.080738,  -0.053092,
                               0.027281, 0.059692, 0.051299,  0.071658,  0.015745,  -0.0099817, 0.026852, 0.034215, -0.037346, -0.025099};
    float coordinates_y[20] = {-0.052672, -0.072924, -0.054435, -0.081549, -0.0025696, -0.071839, -0.049375, 0.062652,  -0.0067283, 0.062291,
                               -0.036252, -0.051579, 0.011029,  0.0039782, 0.059389,   0.038869,  -0.030329, -0.014595, -0.014957,  0.013217};
    int peakCount = 20;
    detectorPeaks_m_t detectorPeaks_m = {coordinates_x, coordinates_y, peakCount};

    const int maxAssambledLatticesCount = 2;
    Lattice_t assembledLattices[maxAssambledLatticesCount];
    int assembledLatticesCount;
    IndexerPlain_index(indexer, assembledLattices, &assembledLatticesCount, maxAssambledLatticesCount, &detectorPeaks_m);

    printf("assembledLatticesCount: %d\n\na: %f %f %f\nb: %f %f %f\nc: %f %f %f\n", assembledLatticesCount, assembledLattices[0].ax, assembledLattices[0].ay,
           assembledLattices[0].az, assembledLattices[0].bx, assembledLattices[0].by, assembledLattices[0].bz, assembledLattices[0].cx, assembledLattices[0].cy,
           assembledLattices[0].cz);


    IndexerPlain_delete(indexer);
    ExperimentSettings_delete(experimentSettings);
}


void test_filterSamplePointsForNorm()
{
    vector<Matrix3Xf> samplePoints(100);
    for (int i = 0; i < 100; ++i)
    {
        loadEigenMatrixFromDisk(samplePoints[i], "workfolder/samplePoints");
    }

    ArrayXf allowedNorms;
    loadEigenMatrixFromDisk(allowedNorms, "workfolder/allowedNorms");

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; ++i)
    {
        filterSamplePointsForNorm(samplePoints[i], allowedNorms, 0.03);
    }
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "duration: " << duration << "ms" << endl << endl;

    ofstream out1("workfolder/samplePointsFiltered");
    out1 << samplePoints[0];
}

void test_indexerAutocorrPrefit()
{
    ExperimentSettings experimentSettings = getExperimentSettingLys();

    IndexerAutocorrPrefit indexer(experimentSettings);

    stringstream ss;
    int runNumber = 0;
    chrono::duration<int64_t, milli>::rep totalDuration(0);
    try
    {
        while (1)
        {
            runNumber++;

            Matrix2Xf detectorPeaks_m;
            ss.str("");
            ss.clear();
            ss << "workfolder/detectorPeaks_m__run" << runNumber;
            loadEigenMatrixFromDisk(detectorPeaks_m, ss.str());

            cout << "runNumber " << runNumber << endl;

            vector<Lattice> assembledLattices;

            chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
            indexer.index(assembledLattices, detectorPeaks_m);
            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
            cout << "duration: " << duration << "ms" << endl << endl;
            totalDuration += duration;

            ss.str("");
            ss.clear();
            ss << "workfolder/lattices__run" << runNumber;
            ofstream outfile(ss.str());
            for (uint16_t i = 0; i < assembledLattices.size(); ++i)
            {
                outfile << assembledLattices[i] << endl << endl;
            }
            outfile.close();
        }
    }
    catch (CustomException& e)
    {
        cout << e.what();

        cout << "custom exception caught" << endl;
        cout << "no more files left" << endl << endl;
        cout << "total duration: " << totalDuration << endl << endl;
    }
    catch (exception& e)
    {
        cout << e.what();

        cout << "general exception caught" << endl;
    }
}

void test_indexerPlain()
{
    ExperimentSettings experimentSettings = getExperimentSettingLys();

    IndexerPlain indexer(experimentSettings);
    indexer.setSamplingPitch(IndexerPlain::SamplingPitch::standardWithSeondaryMillerIndices);

    stringstream ss;
    int runNumber = 0;
    chrono::duration<int64_t, milli>::rep totalDuration(0);
    try
    {
        while (1)
        {
            runNumber++;

            Matrix2Xf detectorPeaks_m;
            ss.str("");
            ss.clear();
            ss << "workfolder/detectorPeaks_m__run" << runNumber;
            loadEigenMatrixFromDisk(detectorPeaks_m, ss.str());

            cout << "runNumber " << runNumber << endl;

            vector<Lattice> assembledLattices;

            chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
            indexer.index(assembledLattices, detectorPeaks_m);
            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
            cout << "duration: " << duration << "ms" << endl << endl;
            totalDuration += duration;

            ss.str("");
            ss.clear();
            ss << "workfolder/lattices__run" << runNumber;
            ofstream outfile(ss.str());
            for (uint16_t i = 0; i < assembledLattices.size(); ++i)
            {
                outfile << assembledLattices[i] << endl << endl;
            }
            outfile.close();
        }
    }
    catch (CustomException& e)
    {
        cout << e.what();

        cout << "custom exception caught" << endl;
        cout << "no more files left" << endl << endl;
        cout << "total duration: " << totalDuration << endl << endl;
    }
    catch (exception& e)
    {
        cout << e.what();

        cout << "general exception caught" << endl;
    }
}

void test_dbscan()
{
    Matrix3Xf points;

    loadEigenMatrixFromDisk(points, "workfolder/autocorrelationPoints");

    float maxEpsilon = 0.0037;
    float maxPossiblePointNorm = 1;
    Dbscan dbscan(maxEpsilon, maxPossiblePointNorm);

    uint16_t minPoints = 2;
    float epsilon = 0.0037;

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    vector<Dbscan::cluster_t> clusters;
    dbscan.computeClusters(clusters, points, minPoints, epsilon);

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "duration: " << duration << "ms" << endl << endl;

    RowVectorXf clusterIndex = RowVectorXf::Zero(points.cols());
    for (uint32_t i = 0; i < clusters.size(); ++i)
    {
        const auto& cluster = clusters[i];
        for (uint32_t index : cluster)
        {
            if (clusterIndex[index] != 0)
            {
                cerr << "node belonging to two clusters!";
            }
            clusterIndex[index] = i + 1;
        }
    }

    cout << "clusters found: " << clusters.size();

    ofstream out1("workfolder/clusterIndex");
    out1 << clusterIndex;
}

void test_pointAutocorrelation()
{
    Matrix3Xf autocorrelationPoints;
    VectorXi centerPointIndices;
    VectorXi shiftedPointIndices;
    Matrix3Xf points;
    float maxNormInAutocorrelation = 0.05;
    float minNormInAutocorrelation = 0.02;

    loadEigenMatrixFromDisk(points, "workfolder/points");

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    //    getPointAutocorrelation(autocorrelationPoints, points, maxNormInAutocorrelation);
    getPointAutocorrelation(autocorrelationPoints, centerPointIndices, shiftedPointIndices, points, minNormInAutocorrelation, maxNormInAutocorrelation);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "duration: " << duration << "ms" << endl << endl;

    ofstream out1("workfolder/autocorrelationPoints");
    out1 << autocorrelationPoints;
    ofstream out2("workfolder/centerPointIndices");
    out2 << centerPointIndices;
    ofstream out3("workfolder/shiftedPointIndices");
    out3 << shiftedPointIndices;
}

void test_latticeAssembler()
{
    LatticeAssembler::accuracyConstants_t accuracyConstants_LatticeAssembler;

    accuracyConstants_LatticeAssembler.maxCountGlobalPassingWeightFilter = 500;
    accuracyConstants_LatticeAssembler.maxCountLocalPassingWeightFilter = 15;
    accuracyConstants_LatticeAssembler.maxCountPassingRelativeDefectFilter = 50;
    accuracyConstants_LatticeAssembler.minPointsOnLattice = 5;

    Vector2f determinantRange(196748.751775786, 295123.127663679);

    LatticeAssembler latticeAssembler(determinantRange, accuracyConstants_LatticeAssembler);

    Matrix3Xf candidateVectors;
    RowVectorXf candidateVectorWeights;
    vector<std::vector<uint16_t>> pointIndicesOnVectors;
    Matrix3Xf pointsToFitInReciprocalSpace;

    loadEigenMatrixFromDisk(candidateVectors, "workfolder/candidateVectors");
    loadEigenMatrixFromDisk(candidateVectorWeights, "workfolder/candidateVectorWeights");
    pointIndicesOnVectors.reserve(candidateVectorWeights.size());
    for (uint16_t i = 0; i < candidateVectorWeights.size(); ++i)
    {
        stringstream pathStream;
        pathStream << "workfolder/pointIndicesOnVector" << i;
        ifstream file(pathStream.str());
        istream_iterator<uint16_t> startFile(file), end;
        pointIndicesOnVectors.emplace_back(startFile, end);
    }
    loadEigenMatrixFromDisk(pointsToFitInReciprocalSpace, "workfolder/pointsToFitInReciprocalSpace");

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    vector<Lattice> assembledLattices;
    vector<LatticeAssembler::assembledLatticeStatistics_t> assembledLatticesStatistics;
    latticeAssembler.assembleLattices(assembledLattices, assembledLatticesStatistics, candidateVectors, candidateVectorWeights, pointIndicesOnVectors,
                                      pointsToFitInReciprocalSpace);

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "duration: " << duration << "ms" << endl << endl;

    for (uint16_t i = 0; i < assembledLattices.size(); ++i)
    {
        cout << assembledLattices[i] << endl << endl;

        cout << "meanDefect: " << assembledLatticesStatistics[i].meanDefect << endl
             << "meanRelativeDefect: " << assembledLatticesStatistics[i].meanRelativeDefect << endl
             << "occupiedLatticePointsCount: " << assembledLatticesStatistics[i].occupiedLatticePointsCount << endl
             << endl;
    }
}

void test_sparsePeakFinder()
{
    Matrix3Xf pointPositions;
    RowVectorXf pointValues;

    loadEigenMatrixFromDisk(pointPositions, "workfolder/pointPositions");
    loadEigenMatrixFromDisk(pointValues, "workfolder/pointValues");

    float minDistanceBetweenRealPeaks = 25;
    float maxPossiblePointNorm = 120;
    SparsePeakFinder sparsePeakFinder(minDistanceBetweenRealPeaks, maxPossiblePointNorm);

    std::ofstream ofs2("workfolder/pointValues_cpp", std::ofstream::out);
    ofs2 << pointValues.transpose().eval();

    std::ofstream ofs3("workfolder/pointPositions_cpp", std::ofstream::out);
    ofs3 << pointPositions.transpose().eval();

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    sparsePeakFinder.findPeaks_fast(pointPositions, pointValues);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "duration: " << duration << "ms" << endl;

    std::ofstream ofs("workfolder/peakPositions", std::ofstream::out);
    ofs << pointPositions.transpose().eval();

    std::ofstream ofs1("workfolder/peakValues", std::ofstream::out);
    ofs1 << pointValues.transpose().eval();
}

void test_hillClimbing()
{
    Matrix3Xf pointsToTransform;
    Matrix3Xf positionsToOptimize;

    float unitPitch = 0.05;
    float tolerance = 0.02;
    VectorXf radii = (VectorXf(2) << 38.2457, 80.2551).finished();
    SamplePointsGenerator samplePointsGenerator;
    samplePointsGenerator.getTightGrid(positionsToOptimize, unitPitch, tolerance, radii);
    loadEigenMatrixFromDisk(pointsToTransform, "workfolder/positionsToTransform");

    //    positionsToOptimize = positionsToOptimize.col(0).eval();

    HillClimbingOptimizer optimizer;

    HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbingOptimizer_accuracyConstants;

    hillClimbingOptimizer_accuracyConstants.functionSelection = 1;
    hillClimbingOptimizer_accuracyConstants.optionalFunctionArgument = 1;
    hillClimbingOptimizer_accuracyConstants.maxCloseToPointDeviation = 0.15;

    hillClimbingOptimizer_accuracyConstants.initialIterationCount = 40;
    hillClimbingOptimizer_accuracyConstants.calmDownIterationCount = 5;
    hillClimbingOptimizer_accuracyConstants.calmDownFactor = 0.8;
    hillClimbingOptimizer_accuracyConstants.localFitIterationCount = 8;
    hillClimbingOptimizer_accuracyConstants.localCalmDownIterationCount = 6;
    hillClimbingOptimizer_accuracyConstants.localCalmDownFactor = 0.8;

    hillClimbingOptimizer_accuracyConstants.stepComputationAccuracyConstants.directionChangeFactor = 2.500000000000000;
    hillClimbingOptimizer_accuracyConstants.stepComputationAccuracyConstants.minStep = 0.331259661674998;
    hillClimbingOptimizer_accuracyConstants.stepComputationAccuracyConstants.maxStep = 3.312596616749981;
    hillClimbingOptimizer_accuracyConstants.stepComputationAccuracyConstants.gamma = 0.650000000000000;
    optimizer.setHillClimbingAccuracyConstants(hillClimbingOptimizer_accuracyConstants);

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    optimizer.performOptimization(pointsToTransform, positionsToOptimize);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "duration: " << duration << "ms" << endl;

    std::ofstream ofs("workfolder/optimizedPoints", std::ofstream::out);
    ofs << positionsToOptimize.transpose().eval();

    std::ofstream ofs2("workfolder/lastInverseTransformEvaluation", std::ofstream::out);
    ofs2 << optimizer.getLastInverseTransformEvaluation().transpose().eval();

    std::ofstream ofs3("workfolder/closeToPointsCount", std::ofstream::out);
    ofs3 << optimizer.getCloseToPointsCount().transpose().eval();
}

void test_computeStep()
{
    float maxCloseToPointDeviation = 0.15;

    InverseSpaceTransform t(maxCloseToPointDeviation);
    t.setFunctionSelection(9);
    t.setOptionalFunctionArgument(3);
    t.clearRadialWeightingFlag();
    t.setLocalTransformFlag();

    Matrix3Xf pointsToTransform(3, 5);
    Matrix3Xf positionsToEvaluate(3, 4);

    pointsToTransform << 3.2, 3.3, 3.2, 4.1, 5.12, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
    positionsToEvaluate << 1.7, 2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    t.setPointsToTransform(pointsToTransform);
    t.performTransform(positionsToEvaluate);

    //    cout << t.getInverseTransformEvaluation() << endl << endl << t.getGradient() << endl << endl << t.getCloseToPointsCount() << endl;

    HillClimbingOptimizer h;
    h.previousStepDirection = (Matrix3Xf(3, 4) << 0.5789, 0.6826, 0.3688, 0.6340, 0.4493, 0.0735, 0.8089, 0.3796, 0.6804, 0.7271, 0.4578, 0.6737).finished();
    h.previousStepLength = (Array<float, 1, Eigen::Dynamic>(1, 4) << 0.1, 3, 2, 4).finished();

    HillClimbingOptimizer::stepComputationAccuracyConstants_t stepComputationAccuracyConstants;
    stepComputationAccuracyConstants.directionChangeFactor = 3;
    stepComputationAccuracyConstants.minStep = 0.5;
    stepComputationAccuracyConstants.maxStep = 55;
    stepComputationAccuracyConstants.gamma = 0.2;
    h.setStepComputationAccuracyConstants(stepComputationAccuracyConstants);

    bool useStepOrthogonalization = true;

    h.computeStep(t.getGradient(), t.getCloseToPointsCount(), t.getInverseTransformEvaluation(), useStepOrthogonalization);

    cout << h.step << endl << endl;
}
void test_InverseSpaceTransform()
{
    float maxCloseToPointDeviation = 0.15;

    InverseSpaceTransform t(maxCloseToPointDeviation);
    t.setFunctionSelection(9);
    t.setOptionalFunctionArgument(3);
    t.clearRadialWeightingFlag();
    t.setLocalTransformFlag();

    Matrix3Xf pointsToTransform(3, 5);
    Matrix3Xf positionsToEvaluate(3, 4);

    pointsToTransform << 3.1, 3.1, 3.1, 4.1, 5.1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
    positionsToEvaluate << 1, 2, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    t.setPointsToTransform(pointsToTransform);

    t.performTransform(positionsToEvaluate);

    cout << t.getInverseTransformEvaluation() << endl << endl << t.getGradient() << endl << endl << t.getCloseToPointsCount() << endl;

    std::vector<std::vector<uint16_t>>& peaksCloseToEvaluationPositions_indices = t.getPointsCloseToEvaluationPositions_indices();
    for (uint32_t i = 0; i < peaksCloseToEvaluationPositions_indices.size(); i++)
    {
        cout << endl << i << ": ";
        for (uint32_t j = 0; j < peaksCloseToEvaluationPositions_indices[i].size(); j++)
        {
            cout << peaksCloseToEvaluationPositions_indices[i][j];
        }
    }
}

static ExperimentSettings getExperimentSettingLys()
{
    float coffset_m = 0.567855;
    float clen_mm = -439.9992;
    float beamEenergy_eV = 8.0010e+03;
    float divergenceAngle_deg = 0.05 * M_PI / 180;
    float nonMonochromaticity = 0.005;
    float pixelLength_m = 110e-6;
    float detectorRadius_pixel = 750;

    Vector3f a_star(+0.0945252, -0.0433391, +0.0644485);
    Vector3f b_star(-0.0298714, -0.1177522, -0.0374347);
    Vector3f c_star(+0.1601091, +0.0156280, -0.2065220);
    Matrix3f basis;
    basis << a_star, b_star, c_star;
    basis = basis / 10; // nm to A
    Lattice sampleReciprocalLattice_1A(basis);
    float tolerance = 0.02;

    return ExperimentSettings(coffset_m, clen_mm, beamEenergy_eV, divergenceAngle_deg, nonMonochromaticity, pixelLength_m, detectorRadius_pixel,
                              sampleReciprocalLattice_1A, tolerance);
}
