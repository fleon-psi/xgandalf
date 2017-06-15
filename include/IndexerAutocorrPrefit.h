/*
 * IndexerAutocorrPrefit.h
 *
 *  Created on: 08.06.2017
 *      Author: Yaro
 */

#ifndef INDEXERAUTOCORRPREFIT_H_
#define INDEXERAUTOCORRPREFIT_H_

#include <IndexerBase.h>

class IndexerAutocorrPrefit: public IndexerBase {
public:
    enum class SamplingPitch {
        extremelyLoose,
        loose,
        standard,
        dense,
        extremelyDense
    };

    IndexerAutocorrPrefit(const ExperimentSettings& experimentSettings);
    IndexerAutocorrPrefit(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath);

    void index(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m);

    void setSamplingPitch(SamplingPitch samplingPitch);
    void setSamplingPitch(float unitPitch);

private:
    void precompute();

    void getGoodAutocorrelationPoints(Eigen::Matrix3Xf& goodAutocorrelationPoints, Eigen::RowVectorXf& goodAutocorrelationPointWeights,
            const Eigen::Matrix3Xf& points, uint32_t maxAutocorrelationPointsCount);
    void autocorrPrefit(const Eigen::Matrix3Xf& reciprocalPeaks_A, Eigen::Matrix3Xf& samplePoints,
            HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_autocorr);

    Eigen::Matrix3Xf precomputedSamplePoints;

    HillClimbingOptimizer hillClimbingOptimizer;
    SparsePeakFinder sparsePeakFinder;
    InverseSpaceTransform inverseSpaceTransform;

    float maxCloseToPeakDeviation;
    float maxNormInAutocorrelation;
    float minNormInAutocorrelation;
    float dbscanEpsilon;
    Dbscan dbscan;
};

#endif /* INDEXERAUTOCORRPREFIT_H_ */
