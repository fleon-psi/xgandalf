#ifndef ADAPTIONS_CRYSTFEL_SIMPLE_DIFFRACTION_PATTERN_PREDICTION_H
#define ADAPTIONS_CRYSTFEL_SIMPLE_DIFFRACTION_PATTERN_PREDICTION_H

#include "ExperimentSettings.h"
#include "indexerData.h"
#include "projectionData.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SimpleDiffractionPatternPrediction SimpleDiffractionPatternPrediction;

SimpleDiffractionPatternPrediction* SimpleDiffractionPrediction_new(ExperimentSettings* experimentSettings);
void SimpleDiffractionPatternPrediction_delete(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction);

void SDPP_getPeaksOnEwaldSphere(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction, reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A,
                           Lattice_t lattice);
void SDPP_predictPattern(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction, millerIndices_t* millerIndices,
                    projectionDirections_t* projectionDirections, Lattice_t lattice);

#ifdef __cplusplus
}
#endif


#endif
