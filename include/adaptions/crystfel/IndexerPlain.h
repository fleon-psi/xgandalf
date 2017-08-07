
#ifndef ADAPTIONS_CRYSTFEL_INDEXER_PLAIN_H
#define ADAPTIONS_CRYSTFEL_INDEXER_PLAIN_H

#include "adaptions/crystfel/ExperimentSettings.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SAMPLING_PITCH_extremelyLoose,
    SAMPLING_PITCH_loose,
    SAMPLING_PITCH_standard,
    SAMPLING_PITCH_dense,
    SAMPLING_PITCH_extremelyDense,

    SAMPLING_PITCH_standardWithSeondaryMillerIndices,
    SAMPLING_PITCH_denseWithSeondaryMillerIndices,
    SAMPLING_PITCH_extremelyDenseWithSeondaryMillerIndices
} samplingPitch_t;

typedef enum GradientDescentIterationsCount {
    GRADIENT_DESCENT_ITERATION_COUNT_verryFew,
    GRADIENT_DESCENT_ITERATION_COUNT_few,
    GRADIENT_DESCENT_ITERATION_COUNT_standard,
    GRADIENT_DESCENT_ITERATION_COUNT_many,
    GRADIENT_DESCENT_ITERATION_COUNT_manyMany,
    GRADIENT_DESCENT_ITERATION_COUNT_extremelyMany,

    GRADIENT_DESCENT_ITERATION_COUNT_custom
} gradientDescentIterationsCount_t;

typedef struct
{
    float* coordinates_x;
    float* coordinates_y;
    int peakCount;
} detectorPeaks_m_t;

typedef struct IndexerPlain IndexerPlain;

IndexerPlain* IndexerPlain_new(ExperimentSettings* experimentSettings);

void IndexerPlain_setSamplingPitch(IndexerPlain* indexerPlain, samplingPitch_t samplingPitch);
void IndexerPlain_setGradientDescentIterationsCount(IndexerPlain* indexerPlain, gradientDescentIterationsCount_t gradientDescentIterationsCount);


void IndexerPlain_index(IndexerPlain* indexerPlain, Lattice_t* assembledLattices, int* assembledLatticesCount, int maxAssambledLatticesCount,
                        const detectorPeaks_m_t* detectorPeaks_m);


#ifdef __cplusplus
}
#endif

#endif