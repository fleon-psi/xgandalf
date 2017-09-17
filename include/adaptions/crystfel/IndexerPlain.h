
#ifndef ADAPTIONS_CRYSTFEL_INDEXER_PLAIN_H
#define ADAPTIONS_CRYSTFEL_INDEXER_PLAIN_H

#include "adaptions/crystfel/ExperimentSettings.h"
#include "adaptions/crystfel/indexerData.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SAMPLING_PITCH_extremelyLoose = 0,
    SAMPLING_PITCH_loose = 1,
    SAMPLING_PITCH_standard = 2,
    SAMPLING_PITCH_dense = 3,
    SAMPLING_PITCH_extremelyDense = 4,

    SAMPLING_PITCH_standardWithSeondaryMillerIndices = 5,
    SAMPLING_PITCH_denseWithSeondaryMillerIndices = 6,
    SAMPLING_PITCH_extremelyDenseWithSeondaryMillerIndices = 7,

    SAMPLING_PITCH_lastEnum
} samplingPitch_t;

typedef enum {
    GRADIENT_DESCENT_ITERATION_COUNT_verryFew = 0,
    GRADIENT_DESCENT_ITERATION_COUNT_few = 1,
    GRADIENT_DESCENT_ITERATION_COUNT_standard = 2,
    GRADIENT_DESCENT_ITERATION_COUNT_many = 3,
    GRADIENT_DESCENT_ITERATION_COUNT_manyMany = 4,
    GRADIENT_DESCENT_ITERATION_COUNT_extremelyMany = 5,

    GRADIENT_DESCENT_ITERATION_COUNT_custom = 6,

    GRADIENT_DESCENT_ITERATION_COUNT_lastEnum
} gradientDescentIterationsCount_t;


typedef struct IndexerPlain IndexerPlain;

IndexerPlain* IndexerPlain_new(ExperimentSettings* experimentSettings);
void IndexerPlain_delete(IndexerPlain* indexerPlain);

void IndexerPlain_setSamplingPitch(IndexerPlain* indexerPlain, samplingPitch_t samplingPitch);
void IndexerPlain_setGradientDescentIterationsCount(IndexerPlain* indexerPlain, gradientDescentIterationsCount_t gradientDescentIterationsCount);
void IndexerPlain_setRefineWithExactLattice(IndexerPlain* indexerPlain, int flag);

void IndexerPlain_index(IndexerPlain* indexerPlain, Lattice_t* assembledLattices, int* assembledLatticesCount, int maxAssambledLatticesCount,
                        reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A);

void backProjectDetectorPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, const ExperimentSettings* experimentSettings, const float* coordinates_x,
                              const float* coordinates_y, int peakCount);

void reorderLattice(const Lattice_t* prototype, Lattice_t* lattice);
void reduceLattice(Lattice_t* lattice);


#ifdef __cplusplus
}
#endif

#endif