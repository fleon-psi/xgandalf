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

    void index_standard(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m);

private:
    void construct();

    void precomputeIndexingStrategy_balanced();
    void precomputeSamplePoints_balanced();

    void clearSamplePointsWithLowInverseFunctionEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation,
            float minFunctionEvaluation);
    void filterSamplePointsForInverseFunctionEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation, uint32_t maxToTakeCount);

    ExperimentSettings experimentSettings;
    SamplePointsGenerator samplePointsGenerator;
    DetectorToReciprocalSpaceTransform detectorToReciprocalSpaceTransform;
    HillClimbingOptimizer hillClimbingOptimizer;

    Eigen::Matrix3Xf samplePoints_standard;
    SparsePeakFinder sparsePeakFinder_standard;
    InverseSpaceTransform inverseSpaceTransform_standard;
    
    LatticeAssembler latticeAssembler;
    
    //just for less reallocation
    std::vector< uint32_t > sortIndices;  //to avoid frequent reallocation

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* INDEXER_H_ */
