//============================================================================
// Name        : indexer.cpp
// Author      : Yaro
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream> 
#include <chrono>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include "Lattice.h"
#include "SamplePointsGenerator.h"
#include "InverseSpaceTransform.h"
#include "HillClimbingOptimizer.h"
#include "eigenDiskImport.h"
#include "SparsePeakFinder.h"
#include "LatticeAssembler.h"
#include "tests.h"
#include "pointAutocorrelation.h"

using namespace std;
using namespace Eigen;

int main()
{
    test_indexerAutocorrPrefit();

    return 0;
}

