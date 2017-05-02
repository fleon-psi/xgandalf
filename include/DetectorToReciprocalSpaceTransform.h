/*
 * DetectorToReciprocalSpaceTransform.h
 *
 *  Created on: 01.05.2017
 *      Author: Yaro
 */

#ifndef DETECTORTORECIPROCALSPACETRANSFORM_H_
#define DETECTORTORECIPROCALSPACETRANSFORM_H_

#include <Eigen/Dense>
#include <ExperimentSettings.h>

class DetectorToReciprocalSpaceTransform {
public:
    DetectorToReciprocalSpaceTransform(const ExperimentSettings& experimentSettings);
        
    //coordinate system same as reciprocal x-z
    void computeReciprocalPeaksFromDetectorPeaks(Eigen::Matrix3Xf& reciprocalPeaks_A, const Eigen::Matrix2Xf& detectorPeaks_m );
    
private:
    float reciprocal_lambda_1A;
    float detectorDistance_m;
};

#endif /* DETECTORTORECIPROCALSPACETRANSFORM_H_ */
