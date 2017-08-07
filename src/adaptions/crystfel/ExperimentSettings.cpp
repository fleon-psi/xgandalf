#include "adaptions/crystfel/ExperimentSettings.h"
#include "ExperimentSettings.h"

#include <Eigen/Dense>


ExperimentSettings* ExperimentSettings_new(float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg, float nonMonochromaticity,
                                           const Lattice_t sampleReciprocalLattice_1A, float tolerance)
{
    const Lattice_t& l = sampleReciprocalLattice_1A;

    Eigen::Matrix3f lattice;
    lattice << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;

    return new ExperimentSettings(detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity, lattice, tolerance);
}