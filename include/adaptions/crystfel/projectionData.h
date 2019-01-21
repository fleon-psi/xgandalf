/* 
 * projectionData.h
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
 
#ifndef ADAPTIONS_CRYSTFEL_PROJECTION_DATA_H
#define ADAPTIONS_CRYSTFEL_PROJECTION_DATA_H

#define MAX_PEAK_COUNT_FOR_PROJECTION 80000

typedef struct
{
	float* coordinates_x;
	float* coordinates_y;
	int peakCount;
} detectorPeaks_m_t;

typedef struct
{
	float* coordinates_x;
	float* coordinates_y;
	float* coordinates_z;
	int peakCount;
} projectionDirections_t;

typedef struct
{
	int* h;
	int* k;
	int* l;
	int peakCount;
} millerIndices_t;

#ifdef __cplusplus
extern "C" {
#endif

	void allocDetectorPeaks(detectorPeaks_m_t* detectorPeaks_m);
	void freeDetectorPeaks(detectorPeaks_m_t detectorPeaks_m);

	void allocProjectionDirections(projectionDirections_t* projectionDirections);
	void freeProjectionDirections(projectionDirections_t projectionDirections);

	void allocMillerIndices(millerIndices_t* millerIndices);
	void freeMillerIndices(millerIndices_t millerIndices);

#ifdef __cplusplus
}
#endif


#endif
