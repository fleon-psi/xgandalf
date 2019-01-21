/* 
 * IndexerData.h
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


#ifndef ADAPTIONS_CRYSTFEL_INDEXER_DATA_H
#define ADAPTIONS_CRYSTFEL_INDEXER_DATA_H

#define MAX_PEAK_COUNT_FOR_INDEXER 20000

typedef struct
{
    float* coordinates_x;
    float* coordinates_y;
    float* coordinates_z;
    int peakCount;
} reciprocalPeaks_1_per_A_t;


#ifdef __cplusplus
extern "C" {
#endif

void allocReciprocalPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A);
void freeReciprocalPeaks(reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A);

#ifdef __cplusplus
}
#endif

#endif