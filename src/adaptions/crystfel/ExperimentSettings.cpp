#include "adaptions/crystfel/ExperimentSettings.h"
#include "ExperimentSettings.h"

#include <Eigen/Dense>

ExperimentSettings* ExperimentSettings_new_nolatt(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                                  float nonMonochromaticity, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A,
                                                  float minRealLatticeDeterminant_A3, float maxRealLatticeDeterminant_A3)
{
    return new ExperimentSettings(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity, minRealLatticeVectorLength_A,
                                  maxRealLatticeVectorLength_A, minRealLatticeDeterminant_A3, maxRealLatticeDeterminant_A3);
}

ExperimentSettings* ExperimentSettings_new(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                           float nonMonochromaticity, const Lattice_t sampleReciprocalLattice_1A, float tolerance)
{
    const Lattice_t& l = sampleReciprocalLattice_1A;

    Eigen::Matrix3f lattice;
    lattice << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;

    return new ExperimentSettings(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity, lattice, tolerance);
}

void ExperimentSettings_delete(ExperimentSettings* experimentSettings)
{
    delete experimentSettings;
}