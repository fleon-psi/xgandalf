/*
 * experimentSettings.h
 *
 *  Created on: 09.04.2017
 *      Author: Yaro
 */

#ifndef EXPERIMENTSETTINGS_H_
#define EXPERIMENTSETTINGS_H_

#include "Lattice.h"

class experimentSettings {
public:
    experimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
            const Lattice& sampleReciprocalLattice, float pixelLength_m, float detectorRadius_pixel);

private:
    Lattice sampleRealLattice;
    Lattice sampleReciprocalLattice;

    Eigen::Vector3f realLatticeVectorNorms;
    Eigen::Vector3f realLatticeVectorAngles;
    Eigen::VectorXf differentRealLatticeVectorNorms;
    float realLatticeDeterminant;
    Eigen::Vector3f reciprocalLatticeVectorNorms;
    Eigen::Vector3f reciprocalLatticeVectorAngles;
    float reciprocalLatticeDeterminant;

    float detectorDistance_m;
    float detectorRadius_m;
    float lambda_A, lambdaShort_A, lambdaLong_A;
    float nonMonochromaticity;
    float divergenceAngle_rad;
    float maxResolutionAngle_rad;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* EXPERIMENTSETTINGS_H_ */
