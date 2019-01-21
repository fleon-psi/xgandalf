/* Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
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

cdef extern from "adaptions/crystfel/Lattice.h":
    ctypedef struct Lattice_t:
        float ax
        float ay
        float az

        float bx
        float by
        float bz

        float cx
        float cy
        float cz


cdef extern from "adaptions/crystfel/ExperimentSettings.h":
    ctypedef struct ExperimentSettings:
        pass

    ExperimentSettings* ExperimentSettings_new_nolatt(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                                  float nonMonochromaticity, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A, float reflectionRadius_1_per_A)

    ExperimentSettings* ExperimentSettings_new(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                               float nonMonochromaticity, const Lattice_t sampleReciprocalLattice_1A, float tolerance, float reflectionRadius_1_per_A)

    void ExperimentSettings_delete(ExperimentSettings* experimentSettings)


cdef extern from "adaptions/crystfel/indexerData.h":
    ctypedef struct reciprocalPeaks_1_per_A_t:
        float* coordinates_x
        float* coordinates_y
        float* coordinates_z
        int peakCount

    void allocReciprocalPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A)
    void freeReciprocalPeaks(reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A)

cdef extern from "adaptions/crystfel/projectionData.h":
    ctypedef struct detectorPeaks_m_t:
        float* coordinates_x;
        float* coordinates_y;
        int peakCount;

    ctypedef struct projectionDirections_t:
        float* coordinates_x
        float* coordinates_y
        float* coordinates_z
        int peakCount

    ctypedef struct millerIndices_t:
        int* h
        int* k
        int* l
        int peakCount

    void allocDetectorPeaks(detectorPeaks_m_t* detectorPeaks_m)
    void freeDetectorPeaks(detectorPeaks_m_t detectorPeaks_m)

    void allocProjectionDirections(projectionDirections_t* projectionDirections)
    void freeProjectionDirections(projectionDirections_t projectionDirections)

    void allocMillerIndices(millerIndices_t* millerIndices)
    void freeMillerIndices(millerIndices_t millerIndices)

cdef extern from "adaptions/crystfel/SimpleDiffractionPatternPrediction.h":
    ctypedef struct SimpleDiffractionPatternPrediction:
        pass

    SimpleDiffractionPatternPrediction* SimpleDiffractionPrediction_new(ExperimentSettings* experimentSettings)
    void SimpleDiffractionPatternPrediction_delete(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction)

    void SDPP_getPeaksOnEwaldSphere(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction, reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A,
                               Lattice_t lattice)
    void SDPP_predictPattern(SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction, millerIndices_t* millerIndices,
                        projectionDirections_t* projectionDirections, Lattice_t lattice)


cdef extern from "adaptions/crystfel/IndexerPlain.h":
    ctypedef enum samplingPitch_t:
        SAMPLING_PITCH_extremelyLoose = 0
        SAMPLING_PITCH_loose = 1
        SAMPLING_PITCH_standard = 2
        SAMPLING_PITCH_dense = 3
        SAMPLING_PITCH_extremelyDense = 4

        SAMPLING_PITCH_standardWithSeondaryMillerIndices = 5
        SAMPLING_PITCH_denseWithSeondaryMillerIndices = 6
        SAMPLING_PITCH_extremelyDenseWithSeondaryMillerIndices = 7

        SAMPLING_PITCH_lastEnum


    ctypedef enum gradientDescentIterationsCount_t:
        GRADIENT_DESCENT_ITERATION_COUNT_verryFew = 0
        GRADIENT_DESCENT_ITERATION_COUNT_few = 1
        GRADIENT_DESCENT_ITERATION_COUNT_standard = 2
        GRADIENT_DESCENT_ITERATION_COUNT_many = 3
        GRADIENT_DESCENT_ITERATION_COUNT_manyMany = 4
        GRADIENT_DESCENT_ITERATION_COUNT_extremelyMany = 5

        GRADIENT_DESCENT_ITERATION_COUNT_custom = 6

        GRADIENT_DESCENT_ITERATION_COUNT_lastEnum

    ctypedef struct IndexerPlain:
        pass


    IndexerPlain* IndexerPlain_new(ExperimentSettings* experimentSettings)
    void IndexerPlain_delete(IndexerPlain* indexerPlain)

    void IndexerPlain_setSamplingPitch(IndexerPlain* indexerPlain, samplingPitch_t samplingPitch)
    void IndexerPlain_setGradientDescentIterationsCount(IndexerPlain* indexerPlain, gradientDescentIterationsCount_t gradientDescentIterationsCount)
    void IndexerPlain_setRefineWithExactLattice(IndexerPlain* indexerPlain, int flag)

    void IndexerPlain_index(IndexerPlain* indexerPlain, Lattice_t* assembledLattices, int* assembledLatticesCount, int maxAssambledLatticesCount,
                            reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A, int* peakCountOnLattices)

    void backProjectDetectorPeaks(reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, const ExperimentSettings* experimentSettings, const float* coordinates_x,
                                  const float* coordinates_y, int peakCount)

    void reorderLattice(const Lattice_t* prototype, Lattice_t* lattice)
    void reduceLattice(Lattice_t* lattice)







