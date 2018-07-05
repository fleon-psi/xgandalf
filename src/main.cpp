//============================================================================
// Name        : indexer.cpp
// Author      : Yaro
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cmath>

#include "HillClimbingOptimizer.h"
#include "InverseSpaceTransform.h"
#include "Lattice.h"
#include "LatticeAssembler.h"
#include "SamplePointsGenerator.h"
#include "SparsePeakFinder.h"
#include "eigenDiskImport.h"
#include "pointAutocorrelation.h"
#include "tests.h"
#include <Eigen/Dense>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Eigen;

int main()
{
    try
    {
		// testPatternPrediction();
        // test_filterSamplePointsForNorm();
        // test_indexerAutocorrPrefit();
        // test_indexerPlain();
        // test_crystfelAdaption();
        // test_crystfelAdaption2();
        // test_latticeReorder();
        // test_gradientDescentRefinement();
        // test_mixedGradientDescentRefinement();
        // test_fixedBasisRefinement();
		// test_fixedBasisRefinementKabsch();
		test_hillClimbing();
    }
    catch (exception& e)
    {
        cout << e.what();
    }

	cout << endl << "done";
    getchar();

    return 0;
}
