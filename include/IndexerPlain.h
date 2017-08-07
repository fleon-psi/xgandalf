/*
 * IndexerPlain.h
 *
 *  Created on: 07.06.2017
 *      Author: Yaro
 */

#ifndef INDEXERPLAIN_H_
#define INDEXERPLAIN_H_

#include "HillClimbingOptimizer.h"
#include <IndexerBase.h>

class IndexerPlain : public IndexerBase
{
  public:
    enum class SamplingPitch
    {
        extremelyLoose,
        loose,
        standard,
        dense,
        extremelyDense,

        standardWithSeondaryMillerIndices,
        denseWithSeondaryMillerIndices,
        extremelyDenseWithSeondaryMillerIndices
    };

    enum class GradientDescentIterationsCount
    {
		exremelyFew,
        few,
        standard,
        many,
        manyMany,
        extremelyMany,

        custom
    };

    IndexerPlain(const ExperimentSettings& experimentSettings);
    IndexerPlain(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath);

    void index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m);

    void setSamplingPitch(SamplingPitch samplingPitch);
    void setSamplingPitch(float unitPitch, bool coverSecondaryMillerIndices);

    void setGradientDescentIterationsCount(GradientDescentIterationsCount gradientDescentIterationsCount);

  private:
    void precompute();

    Eigen::Matrix3Xf precomputedSamplePoints;

    HillClimbingOptimizer hillClimbingOptimizer;
    SparsePeakFinder sparsePeakFinder;
    InverseSpaceTransform inverseSpaceTransform;

    float maxCloseToPointDeviation;

    HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_global;
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_additionalGlobal;
    HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_peaks;
    LatticeAssembler::accuracyConstants_t accuracyConstants_LatticeAssembler;
};

#endif /* INDEXERPLAIN_H_ */
