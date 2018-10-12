#include "adaptions/crystfel/IndexerPlain.h"
#include "IndexerPlain.h"


extern "C" IndexerPlain* IndexerPlain_new(ExperimentSettings* experimentSettings)
{
    return new IndexerPlain(*experimentSettings);
}

extern "C" void IndexerPlain_delete(IndexerPlain* indexerPlain)
{
    delete indexerPlain;
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

void IndexerPlain_setRefineWithExactLattice(IndexerPlain* indexerPlain, int flag)
{
    indexerPlain->setRefineWithExactLattice((bool)flag);
}

extern "C" void IndexerPlain_index(IndexerPlain* indexerPlain, Lattice_t* assembledLattices, int* assembledLatticesCount, int maxAssambledLatticesCount,
                                   reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A, int* peakCountOnLattices)
{
    Eigen::Matrix3Xf reciprocalPeaks_1_per_A_matrix(3, reciprocalPeaks_1_per_A.peakCount);
    for (int i = 0; i < reciprocalPeaks_1_per_A.peakCount; i++)
    {
        reciprocalPeaks_1_per_A_matrix.col(i) << reciprocalPeaks_1_per_A.coordinates_x[i], reciprocalPeaks_1_per_A.coordinates_y[i],
            reciprocalPeaks_1_per_A.coordinates_z[i];
    }

    std::vector<Lattice> assembledLatticesVector;
    std::vector<int> peakCountOnLatticesVector;
    indexerPlain->index(assembledLatticesVector, reciprocalPeaks_1_per_A_matrix, peakCountOnLatticesVector);

    for (*assembledLatticesCount = 0; (int)*assembledLatticesCount < assembledLatticesVector.size() && *assembledLatticesCount < maxAssambledLatticesCount;
         (*assembledLatticesCount)++)
    {
        Eigen::Matrix3f basis = assembledLatticesVector[*assembledLatticesCount].getBasis();
        assembledLattices[*assembledLatticesCount].ax = basis(0, 0);
        assembledLattices[*assembledLatticesCount].ay = basis(1, 0);
        assembledLattices[*assembledLatticesCount].az = basis(2, 0);
        assembledLattices[*assembledLatticesCount].bx = basis(0, 1);
        assembledLattices[*assembledLatticesCount].by = basis(1, 1);
        assembledLattices[*assembledLatticesCount].bz = basis(2, 1);
        assembledLattices[*assembledLatticesCount].cx = basis(0, 2);
        assembledLattices[*assembledLatticesCount].cy = basis(1, 2);
        assembledLattices[*assembledLatticesCount].cz = basis(2, 2);

        if (peakCountOnLattices != NULL)
        {
            peakCountOnLattices[*assembledLatticesCount] = peakCountOnLatticesVector[*assembledLatticesCount];
        }
    }
}


extern "C" void backProjectDetectorPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, const ExperimentSettings* experimentSettings,
                                         const float* coordinates_x, const float* coordinates_y, int peakCount)
{
    DetectorToReciprocalSpaceTransform detectorToReciprocalSpaceTransform(*experimentSettings);
    Eigen::Matrix3Xf reciprocalPeaks_1_per_A_matrix;

    Eigen::Matrix2Xf detectorPeaks_m(2, peakCount);
    for (int i = 0; i < peakCount; i++)
    {
        detectorPeaks_m.col(i) << coordinates_x[i], coordinates_y[i];
    }

    detectorToReciprocalSpaceTransform.computeReciprocalPeaksFromDetectorPeaks(reciprocalPeaks_1_per_A_matrix, detectorPeaks_m);

    int peakNumber;
    for (peakNumber = 0; peakNumber < peakCount && peakNumber < MAX_PEAK_COUNT_FOR_INDEXER; peakNumber++)
    {
        reciprocalPeaks_1_per_A->coordinates_x[peakNumber] = reciprocalPeaks_1_per_A_matrix(0, peakNumber);
        reciprocalPeaks_1_per_A->coordinates_y[peakNumber] = reciprocalPeaks_1_per_A_matrix(1, peakNumber);
        reciprocalPeaks_1_per_A->coordinates_z[peakNumber] = reciprocalPeaks_1_per_A_matrix(2, peakNumber);
    }
    reciprocalPeaks_1_per_A->peakCount = peakNumber;
}

extern "C" void reorderLattice(const Lattice_t* prototype, Lattice_t* lattice)
{
    Eigen::Matrix3f prototypeBasis;
    prototypeBasis << prototype->ax, prototype->bx, prototype->cx, prototype->ay, prototype->by, prototype->cy, prototype->az, prototype->bz, prototype->cz;
    Lattice prototypeLattice(prototypeBasis);

    Eigen::Matrix3f basis;
    basis << lattice->ax, lattice->bx, lattice->cx, lattice->ay, lattice->by, lattice->cy, lattice->az, lattice->bz, lattice->cz;
    Lattice latticeWrapper(basis);

    latticeWrapper.reorder(prototypeLattice.getBasisVectorNorms(), prototypeLattice.getBasisVectorAngles_deg());
    basis = latticeWrapper.getBasis();

    lattice->ax = basis(0, 0);
    lattice->ay = basis(1, 0);
    lattice->az = basis(2, 0);
    lattice->bx = basis(0, 1);
    lattice->by = basis(1, 1);
    lattice->bz = basis(2, 1);
    lattice->cx = basis(0, 2);
    lattice->cy = basis(1, 2);
    lattice->cz = basis(2, 2);
}


void reduceLattice(Lattice_t* lattice, LatticeTransform_t* appliedReductionTransform)
{
    Eigen::Matrix3f basis, reducedBasis;
    basis << lattice->ax, lattice->bx, lattice->cx, lattice->ay, lattice->by, lattice->cy, lattice->az, lattice->bz, lattice->cz;
    Lattice latticeWrapper(basis);

    reducedBasis = latticeWrapper.minimize().getBasis();

    lattice->ax = reducedBasis(0, 0);
    lattice->ay = reducedBasis(1, 0);
    lattice->az = reducedBasis(2, 0);
    lattice->bx = reducedBasis(0, 1);
    lattice->by = reducedBasis(1, 1);
    lattice->bz = reducedBasis(2, 1);
    lattice->cx = reducedBasis(0, 2);
    lattice->cy = reducedBasis(1, 2);
    lattice->cz = reducedBasis(2, 2);

    Eigen::Matrix3f transformMatrix = basis.inverse() * reducedBasis;

    appliedReductionTransform->matrixElement_0_0 = transformMatrix(0, 0);
    appliedReductionTransform->matrixElement_1_0 = transformMatrix(1, 0);
    appliedReductionTransform->matrixElement_2_0 = transformMatrix(2, 0);
    appliedReductionTransform->matrixElement_0_1 = transformMatrix(0, 1);
    appliedReductionTransform->matrixElement_1_1 = transformMatrix(1, 1);
    appliedReductionTransform->matrixElement_2_1 = transformMatrix(2, 1);
    appliedReductionTransform->matrixElement_0_2 = transformMatrix(0, 2);
    appliedReductionTransform->matrixElement_1_2 = transformMatrix(1, 2);
    appliedReductionTransform->matrixElement_2_2 = transformMatrix(2, 2);
}

void restoreLattice(Lattice_t* lattice, LatticeTransform_t* appliedReductionTransform)
{
    Eigen::Matrix3f transformMatrix;

    transformMatrix(0, 0) = appliedReductionTransform->matrixElement_0_0;
    transformMatrix(1, 0) = appliedReductionTransform->matrixElement_1_0;
    transformMatrix(2, 0) = appliedReductionTransform->matrixElement_2_0;
    transformMatrix(0, 1) = appliedReductionTransform->matrixElement_0_1;
    transformMatrix(1, 1) = appliedReductionTransform->matrixElement_1_1;
    transformMatrix(2, 1) = appliedReductionTransform->matrixElement_2_1;
    transformMatrix(0, 2) = appliedReductionTransform->matrixElement_0_2;
    transformMatrix(1, 2) = appliedReductionTransform->matrixElement_1_2;
    transformMatrix(2, 2) = appliedReductionTransform->matrixElement_2_2;

    Eigen::Matrix3f basis, reducedBasis;
    reducedBasis << lattice->ax, lattice->bx, lattice->cx, lattice->ay, lattice->by, lattice->cy, lattice->az, lattice->bz, lattice->cz;

    basis = reducedBasis * transformMatrix.inverse();

    lattice->ax = basis(0, 0);
    lattice->ay = basis(1, 0);
    lattice->az = basis(2, 0);
    lattice->bx = basis(0, 1);
    lattice->by = basis(1, 1);
    lattice->bz = basis(2, 1);
    lattice->cx = basis(0, 2);
    lattice->cy = basis(1, 2);
    lattice->cz = basis(2, 2);
}