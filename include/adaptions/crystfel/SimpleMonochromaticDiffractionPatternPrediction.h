#ifndef ADAPTIONS_CRYSTFEL_SIMPLE_DIFFRACTION_PATTERN_PREDICTION_H
#define ADAPTIONS_CRYSTFEL_SIMPLE_DIFFRACTION_PATTERN_PREDICTION_H

#include "ExperimentSettings.h"
#include "indexerData.h"
#include "projectionData.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SimpleMonochromaticDiffractionPatternPrediction SimpleMonochromaticDiffractionPatternPrediction;

SimpleMonochromaticDiffractionPatternPrediction* SimpleMonochromaticDiffractionPatternPrediction_new(ExperimentSettings* experimentSettings);
void SimpleMonochromaticDiffractionPatternPrediction_delete(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction);

void SMDPP_getPeaksOnEwaldSphere(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction,
                                 reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, Lattice_t lattice);

void SMDPP_predictPattern(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction, millerIndices_t* millerIndices,
                          projectionDirections_t* projectionDirections, Lattice_t lattice);

#ifdef __cplusplus
}
#endif


#endif
