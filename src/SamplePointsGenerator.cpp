/*
 * samplePointsGenerator.cpp
 *
 *  Created on: 10.04.2017
 *      Author: Yaro
 */

#include <SamplePointsGenerator.h>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include "BadInputException.h"
#include "eigenSTLContainers.h"
#include "eigenDiskImport.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

SamplePointsGenerator::SamplePointsGenerator()
{
    precomputedSamplePointsPath = "precomputedSamplePoints";
}
SamplePointsGenerator::SamplePointsGenerator(const std::string& precomputedSamplePointsPath) :
        precomputedSamplePointsPath(precomputedSamplePointsPath)
{
}

inline static float getClosestArrayElement(ArrayXf arr, float value)
{
    int minIndex;
    (arr - value).abs().minCoeff(&minIndex);
    return arr[minIndex];
}

void SamplePointsGenerator::loadPrecomputedSamplePoints(Matrix3Xf& samplePoints, float unitPitch, float tolerance)
{
    Array< float, 1, Eigen::Dynamic > pitches, tolerances;

    stringstream fullPath;
    fullPath << precomputedSamplePointsPath << "/pitches";
    loadEigenMatrixFromDisk(pitches, fullPath.str());

    fullPath.str(string());
    fullPath << precomputedSamplePointsPath << "/tolerances";
    loadEigenMatrixFromDisk(tolerances, fullPath.str());

    fullPath.str(string());
    fullPath << precomputedSamplePointsPath <<
            "/pitch" << getClosestArrayElement(pitches, unitPitch) <<
            "_tolerance" << getClosestArrayElement(tolerances, tolerance);
    MatrixX3f samplePoints_T;
    loadEigenMatrixFromDisk(samplePoints_T, fullPath.str());
    samplePoints = samplePoints_T.transpose();  //workaround, because sample points are stored in a better human readable way

    //old
//    stringstream fullPath;
//    fullPath << precomputedSamplePointsPath << "/pitches";
//    ifstream pitchesFile(fullPath.str());
//    if (!pitchesFile.is_open()) {
//        stringstream errStream;
//        errStream << "Precomputed file " << fullPath.str() << " not found.";
//        throw BadInputException(errStream.str());
//    }
//    istream_iterator< float > startPitchesFile(pitchesFile), end;
//    vector< float > numbers(startPitchesFile, end);
//    pitches = Map< ArrayXf >(&numbers[0], numbers.size());
//
//    fullPath.str(string());
//    numbers.clear();
//    fullPath << precomputedSamplePointsPath << "/tolerances";
//    ifstream tolerancesFile(fullPath.str());
//    if (!tolerancesFile.is_open()) {
//        stringstream errStream;
//        errStream << "Precomputed file " << fullPath.str() << " not found.";
//        throw BadInputException(errStream.str());
//    }
//    istream_iterator< float > startTolerancesFile(tolerancesFile);
//    numbers.assign(startTolerancesFile, end);
//    tolerances = Map< ArrayXf >(&numbers[0], numbers.size());
//
//    fullPath.str(string());
//    numbers.clear();
//    fullPath << precomputedSamplePointsPath <<
//            "/pitch" << getClosestArrayElement(pitches, unitPitch) <<
//            "_tolerance" << getClosestArrayElement(tolerances, tolerance);
//    ifstream samplePointsFile(fullPath.str());
//    if (!samplePointsFile.is_open()) {
//        stringstream errStream;
//        errStream << "Precomputed file " << fullPath.str() << " not found.";
//        throw BadInputException(errStream.str());
//    }
//    istream_iterator< float > startSamplePointsFile(samplePointsFile);
//    numbers.assign(startSamplePointsFile, end);
//    samplePoints = Map< Matrix3Xf >(&numbers[0], 3, numbers.size() / 3); //copy

//    cout << "loading file " << fullPath.str() << endl;
}

void SamplePointsGenerator::getDenseGrid(Matrix3Xf& samplePoints, float unitPitch, float minRadius, float maxRadius)
{
    EigenSTL::vector_Vector3f tmpSamplePoints;

    VectorXf xSamples, ySamples, zSamples;

    int samplesPerRadius = 1 / unitPitch;
    float minRadiusSquared = pow(minRadius, 2);
    float maxRadiusSquared = pow(maxRadius, 2);

    xSamples.setLinSpaced(samplesPerRadius * 2, -maxRadius, maxRadius);
    ySamples.setLinSpaced(samplesPerRadius * 2, -maxRadius, maxRadius);
    zSamples.setLinSpaced(samplesPerRadius, 0, maxRadius);

    for (int xIndex = 0; xIndex < xSamples.size(); xIndex++) {
        for (int yIndex = 0; yIndex < ySamples.size(); yIndex++) {
            for (int zIndex = 0; zIndex < zSamples.size(); zIndex++) {
                Vector3f samplePoint(xSamples[xIndex], ySamples[yIndex], zSamples[zIndex]);
                float squaredNorm = samplePoint.squaredNorm();
                if (squaredNorm >= minRadiusSquared && squaredNorm <= maxRadiusSquared) {
                    tmpSamplePoints.push_back(samplePoint);
                }
            }
        }
    }

    samplePoints = Map< Matrix3Xf >(tmpSamplePoints[0].data(), 3, tmpSamplePoints.size()); //copy
}

void SamplePointsGenerator::getTightGrid(Matrix3Xf& samplePoints, float unitPitch, float tolerance, const VectorXf radii)
{
    for (int i = 0; i < radii.size(); i++) {
        float radius = radii[i];
        float adaptedPitch = unitPitch / radii[i] * radii.maxCoeff();

        Matrix3Xf newSamplingPoints;
        loadPrecomputedSamplePoints(newSamplingPoints, adaptedPitch, tolerance);
        newSamplingPoints *= radius;

        Matrix3Xf tmp(3, samplePoints.cols() + newSamplingPoints.cols());
        tmp << samplePoints, newSamplingPoints;
        samplePoints = tmp;
    }
}
