#include "adaptions/crystfel/IndexerPlain.h"
#include "IndexerPlain.h"


extern "C" IndexerPlain* IndexerPlain_new(ExperimentSettings* experimentSettings)
{
    return new IndexerPlain(*experimentSettings);
}


extern "C" void IndexerPlain_setSamplingPitch(IndexerPlain* indexerPlain, samplingPitch_t samplingPitch)
{
    IndexerPlain::SamplingPitch pitch;
    switch (samplingPitch)
    {
        case SAMPLING_PITCH_extremelyLoose:
            pitch = IndexerPlain::SamplingPitch::extremelyLoose;
            break;
        case SAMPLING_PITCH_loose:
            pitch = IndexerPlain::SamplingPitch::loose;
            break;
        case SAMPLING_PITCH_standard:
            pitch = IndexerPlain::SamplingPitch::standard;
            break;
        case SAMPLING_PITCH_dense:
            pitch = IndexerPlain::SamplingPitch::dense;
            break;
        case SAMPLING_PITCH_extremelyDense:
            pitch = IndexerPlain::SamplingPitch::extremelyDense;
            break;

        case SAMPLING_PITCH_standardWithSeondaryMillerIndices:
            pitch = IndexerPlain::SamplingPitch::standardWithSeondaryMillerIndices;
            break;
        case SAMPLING_PITCH_denseWithSeondaryMillerIndices:
            pitch = IndexerPlain::SamplingPitch::denseWithSeondaryMillerIndices;
            break;
        case SAMPLING_PITCH_extremelyDenseWithSeondaryMillerIndices:
            pitch = IndexerPlain::SamplingPitch::extremelyDenseWithSeondaryMillerIndices;
            break;
        default:
            pitch = IndexerPlain::SamplingPitch::standard;
            break;
    }

    indexerPlain->setSamplingPitch(pitch);
}

extern "C" void IndexerPlain_setGradientDescentIterationsCount(IndexerPlain* indexerPlain, gradientDescentIterationsCount_t gradientDescentIterationsCount)
{
    IndexerPlain::GradientDescentIterationsCount iterationsCount;

    switch (gradientDescentIterationsCount)
    {
        case GRADIENT_DESCENT_ITERATION_COUNT_verryFew:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::exremelyFew;
            break;
        case GRADIENT_DESCENT_ITERATION_COUNT_few:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::few;
            break;
        case GRADIENT_DESCENT_ITERATION_COUNT_standard:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::standard;
            break;
        case GRADIENT_DESCENT_ITERATION_COUNT_many:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::many;
            break;
        case GRADIENT_DESCENT_ITERATION_COUNT_manyMany:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::manyMany;
            break;
        case GRADIENT_DESCENT_ITERATION_COUNT_extremelyMany:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::extremelyMany;
            break;

        default:
            iterationsCount = IndexerPlain::GradientDescentIterationsCount::standard;
            break;
    }

    indexerPlain->setGradientDescentIterationsCount(iterationsCount);
}


extern "C" void IndexerPlain_index(IndexerPlain* indexerPlain, Lattice_t* assembledLattices, int* assembledLatticesCount, int maxAssambledLatticesCount,
           const detectorPeaks_m_t* detectorPeaks_m)
{
    Eigen::Matrix2Xf detectorPeaks_m_matrix(2, detectorPeaks_m->peakCount);
    for (int i = 0; i < detectorPeaks_m->peakCount; i++)
    {
        detectorPeaks_m_matrix.col(i) << detectorPeaks_m->coordinates_x[i], detectorPeaks_m->coordinates_y[i];
    }

    std::vector<Lattice> assembledLatticesVector;
    indexerPlain->index(assembledLatticesVector, detectorPeaks_m_matrix);

    for (*assembledLatticesCount = 0; *assembledLatticesCount < assembledLatticesVector.size() && *assembledLatticesCount < maxAssambledLatticesCount;
         (*assembledLatticesCount)++)
    {
        const Eigen::Matrix3f& basis = assembledLatticesVector[*assembledLatticesCount].getBasis();
        assembledLattices[*assembledLatticesCount].ax = basis(0, 0);
        assembledLattices[*assembledLatticesCount].ay = basis(1, 0);
        assembledLattices[*assembledLatticesCount].az = basis(2, 0);
        assembledLattices[*assembledLatticesCount].bx = basis(0, 1);
        assembledLattices[*assembledLatticesCount].by = basis(1, 1);
        assembledLattices[*assembledLatticesCount].bz = basis(2, 1);
        assembledLattices[*assembledLatticesCount].cx = basis(0, 2);
        assembledLattices[*assembledLatticesCount].cy = basis(1, 2);
        assembledLattices[*assembledLatticesCount].cz = basis(2, 2);
    }
}