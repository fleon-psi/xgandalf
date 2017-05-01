/*
 * ReciprocalSpaceTranslation.h
 *
 *  Created on: 01.05.2017
 *      Author: Yaro
 */

#ifndef RECIPROCALSPACETRANSLATION_H_
#define RECIPROCALSPACETRANSLATION_H_

#include <Eigen/Dense>
#include <ExperimentSettings.h>

class ReciprocalSpaceTranslation {
public:
    ReciprocalSpaceTranslation(const ExperimentSettings& experimentSettings);
        
    //coordinate system same as reciprocal x-z
    void computeReciprocalPeaksFromDetectorPeaks(Eigen::Matrix3Xf& reciprocalPeaks_A, const Eigen::Matrix2Xf& detectorPeaks_m );
    
private:
    float reciprocal_lambda_1A;
    float detectorDistance_m;
};

#endif /* RECIPROCALSPACETRANSLATION_H_ */
