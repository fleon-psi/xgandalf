#ifndef ADAPTIONS_CRYSTFEL_EXPERIMENT_SETTINGS_H
#define ADAPTIONS_CRYSTFEL_EXPERIMENT_SETTINGS_H

#include "adaptions/crystfel/Lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ExperimentSettings ExperimentSettings;

ExperimentSettings* ExperimentSettings_new(float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg, float nonMonochromaticity,
                                           const Lattice_t sampleReciprocalLattice_1A, float tolerance);

#ifdef __cplusplus
}
#endif

#endif