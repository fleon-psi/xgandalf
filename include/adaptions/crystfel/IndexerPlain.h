/* 
 * IndexerPlain.h
 * 
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019      Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
 *
 * This file is part of XGANDALF.
 *
 * XGANDALF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of 
 * the License, or (at your option) any later version.
 *
 * XGANDALF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with XGANDALF.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef ADAPTIONS_CRYSTFEL_INDEXER_PLAIN_H
#define ADAPTIONS_CRYSTFEL_INDEXER_PLAIN_H

#include "ExperimentSettings.h"
#include "indexerData.h"

#ifdef __cplusplus
namespace xgandalf {
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
void IndexerPlain_setMaxPeaksToUseForIndexing(IndexerPlain* indexerPlain, int maxPeaksToUseForIndexing);

void IndexerPlain_index(IndexerPlain* indexerPlain, Lattice_t* assembledLattices, int* assembledLatticesCount, int maxAssambledLatticesCount,
                        reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A, int* peakCountOnLattices);

void backProjectDetectorPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, const ExperimentSettings* experimentSettings, const float* coordinates_x,
                              const float* coordinates_y, int peakCount);




#ifdef __cplusplus
}
}
#endif

#endif