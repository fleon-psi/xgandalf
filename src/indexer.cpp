//============================================================================
// Name        : indexer.cpp
// Author      : Yaro
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream> 
#include <Eigen/Dense>
#include "Lattice.h"
#include "SamplePointsGenerator.h"
#include "InverseSpaceTransform.h"
#include "HillClimbingOptimizer.h"
#include "eigenDiskImport.h"

using namespace std;
using namespace Eigen;

int main()
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

    std::ofstream ofs0("workfolder/positionsToOptimize_init", std::ofstream::out);
    ofs0 << positionsToOptimize;

    HillClimbingOptimizer optimizer;

    int functionSelection = 1;
    float optionalFunctionArgument = 1;
    float maxCloseToPeakDeviation = 0.15;
    optimizer.setInverseSpaceTransformAccuracyConstants(functionSelection, optionalFunctionArgument, maxCloseToPeakDeviation);

    int initialIterationCount = 40;
    int calmDownIterationCount = 5;
    float calmDownFactor = 0.8;
    int localFitIterationCount = 8;
    int localCalmDownIterationCount = 6;
    float localCalmDownFactor = 0.8;
    optimizer.setHillClimbingStrategyAccuracyConstants(initialIterationCount, calmDownIterationCount, calmDownFactor, localFitIterationCount,
            localCalmDownIterationCount, localCalmDownFactor);

    float directionChangeFactor = 2.500000000000000;
    float minStep = 0.331259661674998;
    float maxStep = 3.312596616749981;
    float gamma = 0.650000000000000;
    optimizer.setStepComputationAccuracyConstants(gamma, minStep, maxStep, directionChangeFactor);

    optimizer.performOptimization(pointsToTransform, positionsToOptimize);

    std::ofstream ofs("workfolder/optimizedPoints", std::ofstream::out);
    ofs << positionsToOptimize.transpose().eval();

    std::ofstream ofs2("workfolder/lastInverseTransformEvaluation", std::ofstream::out);
    ofs2 << optimizer.getLastInverseTransformEvaluation().transpose().eval();

    return 0;
}

void test_computeStep()
{
    float maxCloseToPeakDeviation = 0.15;

    InverseSpaceTransform t(maxCloseToPeakDeviation);
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

    //    cout << t.getInverseTransformEvaluation() << endl << endl << t.getGradient() << endl << endl << t.getCloseToPeaksCount() << endl;

    HillClimbingOptimizer h;
    h.previousStepDirection = (Matrix3Xf(3, 4) << 0.5789, 0.6826, 0.3688, 0.6340, 0.4493, 0.0735, 0.8089, 0.3796, 0.6804, 0.7271, 0.4578, 0.6737).finished();
    h.previousStepLength = (Array< float, 1, Eigen::Dynamic >(1, 4) << 0.1, 3, 2, 4).finished();
    float gamma = 0.2;
    float minStep = 0.5;
    float maxStep = 55;
    float directionChangeFactor = 3;
    h.setStepComputationAccuracyConstants(gamma, minStep, maxStep, directionChangeFactor);

    bool useStepOrthogonalization = true;

    h.computeStep(t.getGradient(), t.getCloseToPeaksCount(), t.getInverseTransformEvaluation(), useStepOrthogonalization);

    cout << h.step << endl << endl;
}
void test_InverseSpaceTransform()
{
    float maxCloseToPeakDeviation = 0.15;

    InverseSpaceTransform t(maxCloseToPeakDeviation);
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

    cout << t.getInverseTransformEvaluation() << endl << endl << t.getGradient() << endl << endl << t.getCloseToPeaksCount() << endl;

    std::vector< std::vector< uint16_t > >& peaksCloseToEvaluationPositions_indices = t.getPeaksCloseToEvaluationPositions_indices();
    for (uint32_t i = 0; i < peaksCloseToEvaluationPositions_indices.size(); i++) {
        cout << endl << i << ": ";
        for (uint32_t j = 0; j < peaksCloseToEvaluationPositions_indices[i].size(); j++) {
            cout << peaksCloseToEvaluationPositions_indices[i][j];
        }
    }
}
