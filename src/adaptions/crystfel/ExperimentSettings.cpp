/*
 * ExperimentSettings.cpp
 *
 * SimpleMonochromaticDiffractionPatternPrediction.h
 * 
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
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
 
#include "adaptions/crystfel/ExperimentSettings.h"
#include "ExperimentSettings.h"

#include <Eigen/Dense>

ExperimentSettings* ExperimentSettings_new_nolatt(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                                  float nonMonochromaticity, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A,
                                                  float reflectionRadius_1_per_A)
{
    return new ExperimentSettings(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity, minRealLatticeVectorLength_A,
                                  maxRealLatticeVectorLength_A, reflectionRadius_1_per_A);
}

ExperimentSettings* ExperimentSettings_new(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                           float nonMonochromaticity, const Lattice_t sampleReciprocalLattice_1A, float tolerance,
                                           float reflectionRadius_1_per_A)
{
    const Lattice_t& l = sampleReciprocalLattice_1A;

    Eigen::Matrix3f lattice;
    lattice << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;

    return new ExperimentSettings(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity, lattice, tolerance,
                                  reflectionRadius_1_per_A);
}

void ExperimentSettings_delete(ExperimentSettings* experimentSettings)
{
    delete experimentSettings;
}