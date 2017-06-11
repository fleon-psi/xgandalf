/*
 * IndexerPlain.h
 *
 *  Created on: 07.06.2017
 *      Author: Yaro
 */

#ifndef INDEXERPLAIN_H_
#define INDEXERPLAIN_H_

#include <IndexerBase.h>

class IndexerPlain: public IndexerBase {
public:
    IndexerPlain(const ExperimentSettings& experimentSettings);
    IndexerPlain(const ExperimentSettings& experimentSettings, const std::string& precomputedSamplePointsPath);

    void index(std::vector< Lattice >& assembledLattices, const Eigen::Matrix2Xf& detectorPeaks_m);

private:
    void precompute();

    Eigen::Matrix3Xf precomputedSamplePoints;
    
    HillClimbingOptimizer hillClimbingOptimizer;
    SparsePeakFinder sparsePeakFinder;
    InverseSpaceTransform inverseSpaceTransform;
    
    float maxCloseToPeakDeviation;
};

#endif /* INDEXERPLAIN_H_ */
