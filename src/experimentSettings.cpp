/*
 * experimentSettings.cpp
 *
 *  Created on: 09.04.2017
 *      Author: Yaro
 */

#include <experimentSettings.h>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Eigen;

experimentSettings::experimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
        const Lattice& sampleReciprocalLattice, float pixelLength_m, float detectorRadius_pixel) :
        sampleReciprocalLattice(sampleReciprocalLattice), nonMonochromaticity(nonMonochromaticity)
{
    this->sampleReciprocalLattice.minimize();

    sampleRealLattice = this->sampleReciprocalLattice.getReciprocalLattice().minimize();

    realLatticeVectorNorms = sampleRealLattice.getBasisVectorNorms();
    realLatticeVectorAngles = sampleRealLattice.getBasisVectorAngles();
    realLatticeDeterminant = sampleRealLattice.det();
    reciprocalLatticeVectorNorms = this->sampleReciprocalLattice.getBasisVectorNorms();
    reciprocalLatticeVectorAngles = this->sampleReciprocalLattice.getBasisVectorAngles();
    reciprocalLatticeDeterminant = this->sampleReciprocalLattice.det();

    detectorDistance_m = clen_mm * 1e-3 + coffset_m;
    detectorRadius_m = detectorRadius_pixel * pixelLength_m;
    divergenceAngle_rad = divergenceAngle_deg * M_PI / 180;
    maxResolutionAngle_rad = atan(detectorRadius_m / detectorDistance_m);
    ;

    float h_Plank = 4.135667662e-15; //Planck constant [eV*s]
    float c_light = 299792458; //speed of light [m/s]
    lambda_A = h_Plank * c_light / beamEenergy_eV * 1e10;
    lambdaShort_A = lambda_A * (1 - nonMonochromaticity / 2);
    lambdaLong_A = lambda_A * (1 + nonMonochromaticity / 2);

    //norms are ordered due to minimization of matrix
    assert(reciprocalLatticeVectorNorms[0] <= reciprocalLatticeVectorNorms[1] && reciprocalLatticeVectorNorms[1] <= reciprocalLatticeVectorNorms[2]);
    assert(realLatticeVectorNorms[0] <= realLatticeVectorNorms[1] && realLatticeVectorNorms[1] <= realLatticeVectorNorms[2]);

    float minSimilarityFactor = 0.96;
    if (reciprocalLatticeVectorNorms[0] / reciprocalLatticeVectorNorms[1] > minSimilarityFactor) {
        if (reciprocalLatticeVectorNorms[1] / reciprocalLatticeVectorNorms[2] > minSimilarityFactor) {
            differentRealLatticeVectorNorms.resize(1);
            differentRealLatticeVectorNorms[0] = reciprocalLatticeVectorNorms.mean();
        } else {
            differentRealLatticeVectorNorms.resize(2);
            differentRealLatticeVectorNorms[0] = reciprocalLatticeVectorNorms.head(2).mean();
            differentRealLatticeVectorNorms[1] = reciprocalLatticeVectorNorms[2];
        }
    } else if (reciprocalLatticeVectorNorms[1] / reciprocalLatticeVectorNorms[2] > minSimilarityFactor) {
        differentRealLatticeVectorNorms.resize(2);
        differentRealLatticeVectorNorms[0] = reciprocalLatticeVectorNorms[0];
        differentRealLatticeVectorNorms[1] = reciprocalLatticeVectorNorms.tail(2).mean();
    } else {
        differentRealLatticeVectorNorms = reciprocalLatticeVectorNorms;
    }
}

