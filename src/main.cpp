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
#include "Indexer.h"
#include "tests.h"

using namespace std;
using namespace Eigen;

ExperimentSettings getExperimentSettingLys();

int main()
{
    ExperimentSettings experimentSettings = getExperimentSettingLys();

    Indexer indexer(experimentSettings);

    stringstream ss;
    int runNumber = 0;
    try {
        while (1) {
            runNumber++;

            Matrix2Xf detectorPeaks_m;
            ss.str("");
            ss.clear();
            ss << "workfolder/detectorPeaks_m__run" << runNumber;
            loadEigenMatrixFromDisk(detectorPeaks_m, ss.str());

            vector< Lattice > assembledLattices;

            chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
            indexer.index_standard(assembledLattices, detectorPeaks_m);
            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast< chrono::milliseconds >(t2 - t1).count();
//            cout << "duration: " << duration << "ms" << endl << endl;

            ss.str("");
            ss.clear();
            ss << "workfolder/lattices__run" << runNumber;
            ofstream outfile(ss.str());
            for (uint16_t i = 0; i < assembledLattices.size(); ++i) {
                outfile << assembledLattices[i] << endl << endl;
            }
            outfile.close();

            cout << "runNumber " << runNumber << endl;
        }
    } catch (...) {
        cout << "no more files left";
    }

    return 0;
}

ExperimentSettings getExperimentSettingLys()
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
    basis = basis / 10; //nm to A
    Lattice sampleReciprocalLattice_1A(basis);

    return ExperimentSettings(coffset_m, clen_mm, beamEenergy_eV, divergenceAngle_deg, nonMonochromaticity,
            pixelLength_m, detectorRadius_pixel, sampleReciprocalLattice_1A);
}
