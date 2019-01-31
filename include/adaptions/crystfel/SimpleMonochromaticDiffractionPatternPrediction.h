/*
 * SimpleMonochromaticDiffractionPatternPrediction.h
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

#ifndef ADAPTIONS_CRYSTFEL_SIMPLE_DIFFRACTION_PATTERN_PREDICTION_H
#define ADAPTIONS_CRYSTFEL_SIMPLE_DIFFRACTION_PATTERN_PREDICTION_H

#include "ExperimentSettings.h"
#include "indexerData.h"
#include "projectionData.h"

#ifdef __cplusplus
namespace xgandalf
{
    extern "C" {
#endif

    typedef struct SimpleMonochromaticDiffractionPatternPrediction SimpleMonochromaticDiffractionPatternPrediction;

    SimpleMonochromaticDiffractionPatternPrediction* SimpleMonochromaticDiffractionPatternPrediction_new(ExperimentSettings* experimentSettings);
    void
    SimpleMonochromaticDiffractionPatternPrediction_delete(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction);

    void SMDPP_getPeaksOnEwaldSphere(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction,
                                     reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, Lattice_t lattice);

    void SMDPP_predictPattern(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction, millerIndices_t* millerIndices,
                              projectionDirections_t* projectionDirections, Lattice_t lattice);

#ifdef __cplusplus
    }
}
#endif


#endif
