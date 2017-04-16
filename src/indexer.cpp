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
    ArrayXXf x(2, 2);

    float maxCloseToPeakDeviation = 0.15;

    InverseSpaceTransform t(maxCloseToPeakDeviation);
    t.setFunctionSelection(9);
    t.setOptionalFunctionArgument(1);
    t.setRadialWeightingFlag();
    t.clearLocalTransformFlag();

    Matrix3Xf pointsToTransform(3, 5);
    Matrix3Xf positionsToEvaluate(3, 4);

    pointsToTransform << 1.1, 2.1, 3.1, 4.1, 5.1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
    positionsToEvaluate << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    t.performTransform(pointsToTransform, positionsToEvaluate);

    cout << t.getInverseTransformEvaluation() << endl << endl << t.getGradient() << endl << endl << t.getCloseToPeaksCount();

    return 0;
}
