//============================================================================
// Name        : indexer.cpp
// Author      : Yaro
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <Eigen/Dense>
#include "Lattice.h"
#include "SamplePointsGenerator.h"
#include "InverseSpaceTransform.h"

using namespace std;
using namespace Eigen;

int main()
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

    return 0;
}

void testInverseSpaceTransform()
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
