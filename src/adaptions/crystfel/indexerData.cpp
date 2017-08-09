#include "adaptions/crystfel/indexerData.h"


void allocReciprocalPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A)
{
    reciprocalPeaks_1_per_A->coordinates_x = new float[MAX_PEAK_COUNT_FOR_INDEXER];
    reciprocalPeaks_1_per_A->coordinates_y = new float[MAX_PEAK_COUNT_FOR_INDEXER];
    reciprocalPeaks_1_per_A->coordinates_z = new float[MAX_PEAK_COUNT_FOR_INDEXER];
}

void freeReciprocalPeaks(reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A)
{
    delete[] reciprocalPeaks_1_per_A.coordinates_x;
    delete[] reciprocalPeaks_1_per_A.coordinates_y;
    delete[] reciprocalPeaks_1_per_A.coordinates_z;
}