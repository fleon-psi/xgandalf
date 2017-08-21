#ifndef ADAPTIONS_CRYSTFEL_EXPERIMENT_SETTINGS_H
#define ADAPTIONS_CRYSTFEL_EXPERIMENT_SETTINGS_H

#include "adaptions/crystfel/Lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ExperimentSettings ExperimentSettings;

ExperimentSettings* ExperimentSettings_new_nolatt(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                                  float nonMonochromaticity, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A);

ExperimentSettings* ExperimentSettings_new(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                           float nonMonochromaticity, const Lattice_t sampleReciprocalLattice_1A, float tolerance);

void ExperimentSettings_delete(ExperimentSettings* experimentSettings);

#ifdef __cplusplus
}
#endif

#endif