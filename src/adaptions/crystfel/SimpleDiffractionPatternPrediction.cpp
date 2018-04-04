#include "adaptions/crystfel/SimpleDiffractionPatternPrediction.h"
#include "SimpleDiffractionPatternPrediction.h"

extern "C" SimpleDiffractionPatternPrediction* SimpleDiffractionPrediction_new(ExperimentSettings* experimentSettings)
{
    return new SimpleDiffractionPatternPrediction(*experimentSettings);
}

extern "C" void SimpleDiffractionPatternPrediction_delete(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction)
{
    delete simpleDiffractionPatternPrediction;
}

extern "C" void SDPP_getPeaksOnEwaldSphere(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction,
                                      reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, Lattice_t lattice)
{
    Eigen::Matrix3Xf reciprocalPeaks_1_per_A_matrix;
    Eigen::Matrix3Xi millerIndices_matrix;

    const Lattice_t& l = lattice;
    Eigen::Matrix3f basis;
    basis << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;
    Lattice lattice_class(basis);

    simpleDiffractionPatternPrediction->getPeaksOnEwaldSphere(reciprocalPeaks_1_per_A_matrix, millerIndices_matrix, lattice_class);

    int peakCount = std::min(MAX_PEAK_COUNT_FOR_PROJECTION, (int)reciprocalPeaks_1_per_A_matrix.cols());
    reciprocalPeaks_1_per_A->peakCount = peakCount;
    for (int i = 0; i < peakCount; i++)
    {
        reciprocalPeaks_1_per_A->coordinates_x[i] = reciprocalPeaks_1_per_A_matrix(0, i);
        reciprocalPeaks_1_per_A->coordinates_y[i] = reciprocalPeaks_1_per_A_matrix(1, i);
        reciprocalPeaks_1_per_A->coordinates_z[i] = reciprocalPeaks_1_per_A_matrix(2, i);
    }
}

extern "C" void SDPP_predictPattern(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction, millerIndices_t* millerIndices,
                               projectionDirections_t* projectionDirections, Lattice_t lattice)
{
    Eigen::Matrix2Xf predictedPeaks_m_matrix;
    Eigen::Matrix3Xi millerIndices_matrix;
    Eigen::Matrix3Xf projectionDirections_matrix;

    const Lattice_t& l = lattice;
    Eigen::Matrix3f basis;
    basis << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;
    Lattice lattice_class(basis);

    simpleDiffractionPatternPrediction->predictPattern(predictedPeaks_m_matrix, millerIndices_matrix, projectionDirections_matrix, lattice_class);

    int peakCount = std::min(MAX_PEAK_COUNT_FOR_PROJECTION, (int)predictedPeaks_m_matrix.cols());
    millerIndices->peakCount = peakCount;
    projectionDirections->peakCount = peakCount;
    for (int i = 0; i < peakCount; i++)
    {
        millerIndices->h[i] = millerIndices_matrix(0, i);
        millerIndices->k[i] = millerIndices_matrix(1, i);
        millerIndices->l[i] = millerIndices_matrix(2, i);

        projectionDirections->coordinates_x[i] = projectionDirections_matrix(0, i);
        projectionDirections->coordinates_y[i] = projectionDirections_matrix(1, i);
        projectionDirections->coordinates_z[i] = projectionDirections_matrix(2, i);
    }
}
