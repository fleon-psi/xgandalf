
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
