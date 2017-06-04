/*
 * Indexer.h
 *
 *  Created on: 30.04.2017
 *      Author: Yaro
 */

#ifndef INDEXER_H_
#define INDEXER_H_

#include <Eigen/Dense>
#include "Lattice.h"
#include "LatticeAssembler.h"
#include "ExperimentSettings.h"
#include "SamplePointsGenerator.h"
#include "DetectorToReciprocalSpaceTransform.h"
#include "HillClimbingOptimizer.h"
#include "SparsePeakFinder.h"

class Indexer {
public:
    Indexer(const ExperimentSettings& experimentSettings);
    Indexer(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath);

    void index_balanced(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m);
    void index_autocorrPrefit(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m);

private:
    void construct();

    void precomputeIndexingStrategy_balanced();
    void precomputeIndexingStrategy_autocorrPrefit();

    void clearSamplePointsWithLowInverseFunctionEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation,
            float minFunctionEvaluation);
    void filterSamplePointsForInverseFunctionEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation, uint32_t maxToTakeCount);

    void autocorrPrefit(Eigen::Matrix3Xf& reciprocalPeaks_A, Eigen::Matrix3Xf& samplePoints);

    ExperimentSettings experimentSettings;
    SamplePointsGenerator samplePointsGenerator;
    DetectorToReciprocalSpaceTransform detectorToReciprocalSpaceTransform;
    HillClimbingOptimizer hillClimbingOptimizer;

    Eigen::Matrix3Xf samplePoints_balanced;
    SparsePeakFinder sparsePeakFinder_balanced;
    Eigen::Matrix3Xf samplePoints_autocorrPrefit;
    SparsePeakFinder sparsePeakFinder_autocorrPrefit;

    LatticeAssembler latticeAssembler;

    //just for less reallocation
    std::vector< uint32_t > sortIndices;  //to avoid frequent reallocation
};

#endif /* INDEXER_H_ */
