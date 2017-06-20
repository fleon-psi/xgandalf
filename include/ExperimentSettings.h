/*
 * experimentSettings.h
 *
 *  Created on: 09.04.2017
 *      Author: Yaro
 */

#ifndef EXPERIMENTSETTINGS_H_
#define EXPERIMENTSETTINGS_H_

#include "Lattice.h"
#include "WrongUsageException.h"

class ExperimentSettings {
public:
    ExperimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
            float pixelLength_m, float detectorRadius_pixel, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A);
    ExperimentSettings(float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg, float nonMonochromaticity,
            float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A);

    ExperimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
            float pixelLength_m, float detectorRadius_pixel, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A,
            float minRealLatticeDeterminant_A3, float maxRealLatticeDeterminant_A3);
    ExperimentSettings(float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg, float nonMonochromaticity,
            float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A, float minRealLatticeDeterminant_A3, float maxRealLatticeDeterminant_A3);

    ExperimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
            float pixelLength_m, float detectorRadius_pixel, const Lattice& sampleReciprocalLattice_1A, float tolerance);
    ExperimentSettings(float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg, float nonMonochromaticity,
            const Lattice& sampleReciprocalLattice_1A, float tolerance);

    float getDetectorDistance_m() const;
    float getDetectorRadius_m() const;
    float getLambda_A() const;
    float getLambdaLong_A() const;
    float getLambdaShort_A() const;
    float getReciprocalLambda_1A() const;
    float getReciprocalLambdaLong_1A() const;
    float getReciprocalLambdaShort_1A() const;
    float getNonMonochromaticity() const;
    float getDivergenceAngle_rad() const;
    float getMaxResolutionAngle_rad() const;

    bool isLatticeParametersKnown() const;

    const Lattice& getSampleRealLattice_A() const;
    const Lattice& getSampleReciprocalLattice_1A() const;
    const Eigen::Array3f& getRealLatticeVectorLengths_A() const;
    const Eigen::Array3f& getRealLatticeVectorAngles_deg() const;
    const Eigen::Array3f& getRealLatticeVectorAnglesNormalized_deg() const;
    float getRealLatticeDeterminant_A3() const;
    const Eigen::Array3f& getReciprocalLatticeVectorLengths_1A() const;
    const Eigen::Array3f& getReciprocalLatticeVectorAngles_deg() const;
    const Eigen::Array3f& getReciprocalLatticeVectorAnglesNormalized_deg() const;
    float getReciprocalLatticeDeterminant_1A3() const;
    float getTolerance() const;

    float getMaxRealLatticeDeterminant_A3() const;
    float getMaxRealLatticeVectorLength_A() const;
    float getMaxReciprocalLatticeDeterminant_1A3() const;
    float getMaxReciprocalLatticeVectorLength_1A() const;
    float getMinRealLatticeDeterminant_A3() const;
    float getMinRealLatticeVectorLength_A() const;
    float getMinReciprocalLatticeDeterminant_1A3() const;
    float getMinReciprocalLatticeVectorLength_1A() const;

    const Eigen::ArrayXf& getDifferentRealLatticeVectorLengths_A() const;

private:
    void constructFromGeometryFileValues(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
            float pixelLength_m, float detectorRadius_pixel);
    void constructFromPrecomputedValues(float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
            float nonMonochromaticity);
    void deduceValuesFromSampleReciprocalLattice();

    float detectorDistance_m;
    float detectorRadius_m;
    float lambda_A, lambdaShort_A, lambdaLong_A;
    float reciprocal_lambda_1A, reciprocal_lambdaShort_1A, reciprocal_lambdaLong_1A;
    float nonMonochromaticity;
    float divergenceAngle_rad;
    float maxResolutionAngle_rad;

    bool latticeParametersKnown;

    float minRealLatticeVectorLength_A;
    float maxRealLatticeVectorLength_A;
    float minRealLatticeDeterminant_A3;
    float maxRealLatticeDeterminant_A3;
    float minReciprocalLatticeVectorLength_1A;
    float maxReciprocalLatticeVectorLength_1A;
    float minReciprocalLatticeDeterminant_1A3;
    float maxReciprocalLatticeDeterminant_1A3;

    //if latticeParametersKnown
    Lattice sampleRealLattice_A;
    Lattice sampleReciprocalLattice_1A;
    Eigen::Array3f realLatticeVectorLengths_A;
    Eigen::Array3f realLatticeVectorAngles_deg;
    Eigen::Array3f realLatticeVectorAnglesNormalized_deg;
    float realLatticeDeterminant_A3;
    Eigen::Array3f reciprocalLatticeVectorLengths_1A;
    Eigen::Array3f reciprocalLatticeVectorAngles_deg;
    Eigen::Array3f reciprocalLatticeVectorAnglesNormalized_deg;
    float reciprocalLatticeDeterminant_1A3;
    float latticeParametersTolerance;

    //if latticeParametersKnown, trivial. if not, set to min and max vector length 
    Eigen::ArrayXf differentRealLatticeVectorLengths_A;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        ;
};

#endif /* EXPERIMENTSETTINGS_H_ */
