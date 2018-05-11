import numpy as np
cimport xgandalf_cpp as cpp


cdef class Xgandalf:
    cdef cpp.ExperimentSettings* experimentSettings

    cdef cpp.IndexerPlain* indexer

    cdef cpp.reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A

    cdef cpp.SimpleDiffractionPatternPrediction* simpleDiffractionPatternPrediction
    cdef cpp.detectorPeaks_m_t detectorPeaks_m
    cdef cpp.projectionDirections_t projectionDirections
    cdef cpp.millerIndices_t millerIndices

    cdef float detectorDistance_m

    def __cinit__(self):
        cpp.allocReciprocalPeaks(&(self.reciprocalPeaks_1_per_A));

        cpp.allocDetectorPeaks(&self.detectorPeaks_m)
        cpp.allocProjectionDirections(&self.projectionDirections)
        cpp.allocMillerIndices(&self.millerIndices)

    def precomputeWithoutLattice(self,
                                    float beamEenergy_eV,
                                    float detectorDistance_m,
                                    float minRealLatticeVectorLength_A,
                                    float maxRealLatticeVectorLength_A,
                                    float tolerance,
                                    int samplingPitch_selector,
                                    int gradientDescentIterationsCount_selector):

        cdef float detectorRadius_m__dummy = 0.5
        cdef float divergenceAngle_deg__dummy = 0.05
        cdef float nonMonochromaticity__dummy = 0.005
        cdef float reflectionRadius_1_per_A__dummy = 0.0002

        self.detectorDistance_m = detectorDistance_m

        self.experimentSettings = cpp.ExperimentSettings_new_nolatt(beamEenergy_eV, detectorDistance_m, detectorRadius_m__dummy, divergenceAngle_deg__dummy,
                                                   nonMonochromaticity__dummy, minRealLatticeVectorLength_A, maxRealLatticeVectorLength_A, reflectionRadius_1_per_A__dummy)

        self.indexer = cpp.IndexerPlain_new(self.experimentSettings)


        cdef cpp.samplingPitch_t samplingPitch = <cpp.samplingPitch_t>samplingPitch_selector
        cdef cpp.gradientDescentIterationsCount_t gradientDescentIterationsCount = <cpp.gradientDescentIterationsCount_t> gradientDescentIterationsCount_selector

        cpp.IndexerPlain_setSamplingPitch(self.indexer, samplingPitch)
        cpp.IndexerPlain_setGradientDescentIterationsCount(self.indexer, gradientDescentIterationsCount)

        self.simpleDiffractionPatternPrediction = cpp.SimpleDiffractionPrediction_new(self.experimentSettings)

    #needs aStar,bStar,cStar be the basis vectors of the primitive reciprocal lattice
    def precomputeWithLattice(self,
                                    float beamEenergy_eV,
                                    float detectorDistance_m,
                                    float[:] aStar,
                                    float[:] bStar,
                                    float[:] cStar,
                                    int useExactLatticeFlag,
                                    float tolerance,
                                    int samplingPitch_selector,
                                    int gradientDescentIterationsCount_selector):

        self.detectorDistance_m = detectorDistance_m

        basis = np.column_stack((aStar,bStar,cStar))
        if( np.linalg.det(basis) < 0 ):
            print("warning!!! Known lattice is not righthanded!!")

        cdef float detectorRadius_m__dummy = 0.09
        cdef float divergenceAngle_deg__dummy = 0.01
        cdef float nonMonochromaticity__dummy = 0.01
        cdef float reflectionRadius_1_per_A__dummy = np.mean(np.linalg.norm(basis, axis=0))*0.05

        cdef cpp.Lattice_t sampleReciprocalLattice_1A
        sampleReciprocalLattice_1A.ax = aStar[0]
        sampleReciprocalLattice_1A.ay = aStar[1]
        sampleReciprocalLattice_1A.az = aStar[2]
        sampleReciprocalLattice_1A.bx = bStar[0]
        sampleReciprocalLattice_1A.by = bStar[1]
        sampleReciprocalLattice_1A.bz = bStar[2]
        sampleReciprocalLattice_1A.cx = cStar[0]
        sampleReciprocalLattice_1A.cy = cStar[1]
        sampleReciprocalLattice_1A.cz = cStar[2]


        self.experimentSettings = cpp.ExperimentSettings_new(beamEenergy_eV, detectorDistance_m, detectorRadius_m__dummy, divergenceAngle_deg__dummy,
                                               nonMonochromaticity__dummy, sampleReciprocalLattice_1A, tolerance, reflectionRadius_1_per_A__dummy)

        self.indexer = cpp.IndexerPlain_new(self.experimentSettings)


        cdef cpp.samplingPitch_t samplingPitch = <cpp.samplingPitch_t>samplingPitch_selector
        cdef cpp.gradientDescentIterationsCount_t gradientDescentIterationsCount = <cpp.gradientDescentIterationsCount_t> gradientDescentIterationsCount_selector

        cpp.IndexerPlain_setSamplingPitch(self.indexer, samplingPitch)
        cpp.IndexerPlain_setGradientDescentIterationsCount(self.indexer, gradientDescentIterationsCount)
        cpp.IndexerPlain_setRefineWithExactLattice(self.indexer, useExactLatticeFlag)

        self.simpleDiffractionPatternPrediction = cpp.SimpleDiffractionPrediction_new(self.experimentSettings)

    def findLattice(self, float[::1] coordinates_x, float[::1] coordinates_y):
        cdef int peakCount = coordinates_x.size
        cpp.backProjectDetectorPeaks(&(self.reciprocalPeaks_1_per_A), self.experimentSettings, &(coordinates_x[0]), &(coordinates_y[0]), peakCount);

        cdef int maxAssambledLatticesCount = 1
        cdef cpp.Lattice_t assembledLattices[1]
        cdef int peakCountOnLattices[1]

        cdef int assembledLatticesCount = 666

        cpp.IndexerPlain_index(self.indexer, &(assembledLattices[0]), &assembledLatticesCount, maxAssambledLatticesCount, self.reciprocalPeaks_1_per_A, &(peakCountOnLattices[0]));

        if assembledLatticesCount > 0:
            a = np.array([assembledLattices[0].ax, assembledLattices[0].ay, assembledLattices[0].az])
            b = np.array([assembledLattices[0].bx, assembledLattices[0].by, assembledLattices[0].bz])
            c = np.array([assembledLattices[0].cx, assembledLattices[0].cy, assembledLattices[0].cz])
            peakCountOnLattice = peakCountOnLattices[0]

            basis = np.column_stack((a,b,c))
            if( np.linalg.det(basis) < 0 ):
                a = -a
                b = -b
                c = -c
        else:
            a = np.array([0,0,0])
            b = np.array([0,0,0])
            c = np.array([0,0,0])
            peakCountOnLattice = 0

        return (a,b,c, peakCountOnLattice)

    def predictPattern(self, aStar, bStar, cStar):
        cdef cpp.Lattice_t lattice_1A
        lattice_1A.ax = aStar[0]
        lattice_1A.ay = aStar[1]
        lattice_1A.az = aStar[2]
        lattice_1A.bx = bStar[0]
        lattice_1A.by = bStar[1]
        lattice_1A.bz = bStar[2]
        lattice_1A.cx = cStar[0]
        lattice_1A.cy = cStar[1]
        lattice_1A.cz = cStar[2]

        #cpp.SDPP_getPeaksOnEwaldSphere(self.simpleDiffractionPatternPrediction, &self.reciprocalPeaks_1_per_A, lattice_1A)
        cpp.SDPP_predictPattern(self.simpleDiffractionPatternPrediction, &self.millerIndices, &self.projectionDirections, lattice_1A)

        peakCount = self.millerIndices.peakCount

        peakCount = self.projectionDirections.peakCount
        cdef float[:] directions_x_view = <float[:peakCount]> self.projectionDirections.coordinates_x
        directions_x = np.asarray(directions_x_view)
        cdef float[:] directions_y_view = <float[:peakCount]> self.projectionDirections.coordinates_y
        directions_y = np.asarray(directions_y_view)
        cdef float[:] directions_z_view = <float[:peakCount]> self.projectionDirections.coordinates_z
        directions_z = np.asarray(directions_z_view)

        detectorCoordinates_x = directions_y/directions_x*self.detectorDistance_m
        detectorCoordinates_y = directions_z/directions_x*self.detectorDistance_m

        return (-detectorCoordinates_x, detectorCoordinates_y)

    def __dealloc__(self):
        cpp.SimpleDiffractionPatternPrediction_delete(self.simpleDiffractionPatternPrediction)
        cpp.ExperimentSettings_delete(self.experimentSettings)
        cpp.IndexerPlain_delete(self.indexer)

        cpp.freeDetectorPeaks(self.detectorPeaks_m)
        cpp.freeProjectionDirections(self.projectionDirections)
        cpp.freeMillerIndices(self.millerIndices)
