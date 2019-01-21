/*
 * projectionData.cpp
 *
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
 
#include "adaptions/crystfel/projectionData.h"

void allocDetectorPeaks(detectorPeaks_m_t* detectorPeaks_m)
{
    detectorPeaks_m->coordinates_x = new float[MAX_PEAK_COUNT_FOR_PROJECTION];
    detectorPeaks_m->coordinates_y = new float[MAX_PEAK_COUNT_FOR_PROJECTION];
}

void freeDetectorPeaks(detectorPeaks_m_t detectorPeaks_m)
{
    delete[] detectorPeaks_m.coordinates_x;
    delete[] detectorPeaks_m.coordinates_y;
}

void allocProjectionDirections(projectionDirections_t* projectionDirections)
{
    projectionDirections->coordinates_x = new float[MAX_PEAK_COUNT_FOR_PROJECTION];
    projectionDirections->coordinates_y = new float[MAX_PEAK_COUNT_FOR_PROJECTION];
    projectionDirections->coordinates_z = new float[MAX_PEAK_COUNT_FOR_PROJECTION];
}

void freeProjectionDirections(projectionDirections_t projectionDirections)
{
    delete[] projectionDirections.coordinates_x;
    delete[] projectionDirections.coordinates_y;
    delete[] projectionDirections.coordinates_z;
}


void allocMillerIndices(millerIndices_t* millerIndices)
{
    millerIndices->h = new int[MAX_PEAK_COUNT_FOR_PROJECTION];
    millerIndices->k = new int[MAX_PEAK_COUNT_FOR_PROJECTION];
    millerIndices->l = new int[MAX_PEAK_COUNT_FOR_PROJECTION];
}

void freeMillerIndices(millerIndices_t millerIndices)
{
    delete[] millerIndices.h;
    delete[] millerIndices.k;
    delete[] millerIndices.l;
}