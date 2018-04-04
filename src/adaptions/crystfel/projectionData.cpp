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