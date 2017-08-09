
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