/*
 * IndexerBase.h
 *
 *  Created on: 06.06.2017
 *      Author: Yaro
 */

#ifndef INDEXERBASE_H_
#define INDEXERBASE_H_

#include "Dbscan.h"
#include "DetectorToReciprocalSpaceTransform.h"
#include "ExperimentSettings.h"
#include "Lattice.h"
#include "LatticeAssembler.h"
#include "SamplePointsGenerator.h"
#include "SparsePeakFinder.h"
#include <Eigen/Dense>

class IndexerBase
{
  public:
    IndexerBase(const ExperimentSettings& experimentSettings);
    IndexerBase(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath);

    virtual ~IndexerBase() = default;

    virtual void index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A) = 0;

  protected:
    void keepSamplePointsWithHighEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation, float minEvaluation);
    void keepSamplePointsWithHighestEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation,
                                               uint32_t maxToTakeCount); // output is sorted

    ExperimentSettings experimentSettings;
    SamplePointsGenerator samplePointsGenerator;

    LatticeAssembler latticeAssembler;

  private:
    std::vector<uint32_t> sortIndices; // to avoid frequent reallocation
};

#endif /* INDEXERBASE_H_ */
