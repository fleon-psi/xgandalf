/*
 * IndexerAutocorrPrefit.cpp
 *
 * SimpleMonochromaticDiffractionPatternPrediction.cpp
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

#include "adaptions/crystfel/SimpleMonochromaticDiffractionPatternPrediction.h"
#include "SimpleMonochromaticDiffractionPatternPrediction.h"

namespace xgandalf
{

    extern "C" SimpleMonochromaticDiffractionPatternPrediction* SimpleMonochromaticDiffractionPatternPrediction_new(ExperimentSettings* experimentSettings)
    {
        return new SimpleMonochromaticDiffractionPatternPrediction(*experimentSettings);
    }

    extern "C" void
    SimpleMonochromaticDiffractionPatternPrediction_delete(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction)
    {
        delete simpleMonochromaticDiffractionPatternPrediction;
    }

    extern "C" void SMDPP_getPeaksOnEwaldSphere(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction,
                                                reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A, Lattice_t lattice)
    {
        Eigen::Matrix3Xf reciprocalPeaks_1_per_A_matrix;
        Eigen::Matrix3Xi millerIndices_matrix;

        const Lattice_t& l = lattice;
        Eigen::Matrix3f basis;
        basis << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;
        Lattice lattice_class(basis);

        simpleMonochromaticDiffractionPatternPrediction->getPeaksOnEwaldSphere(reciprocalPeaks_1_per_A_matrix, millerIndices_matrix, lattice_class);

        int peakCount = std::min(MAX_PEAK_COUNT_FOR_PROJECTION, (int)reciprocalPeaks_1_per_A_matrix.cols());
        reciprocalPeaks_1_per_A->peakCount = peakCount;
        for (int i = 0; i < peakCount; i++)
        {
            reciprocalPeaks_1_per_A->coordinates_x[i] = reciprocalPeaks_1_per_A_matrix(0, i);
            reciprocalPeaks_1_per_A->coordinates_y[i] = reciprocalPeaks_1_per_A_matrix(1, i);
            reciprocalPeaks_1_per_A->coordinates_z[i] = reciprocalPeaks_1_per_A_matrix(2, i);
        }
    }

    extern "C" void SMDPP_predictPattern(SimpleMonochromaticDiffractionPatternPrediction* simpleMonochromaticDiffractionPatternPrediction,
                                         millerIndices_t* millerIndices, projectionDirections_t* projectionDirections, Lattice_t lattice)
    {
        Eigen::Matrix2Xf predictedPeaks_m_matrix;
        Eigen::Matrix3Xi millerIndices_matrix;
        Eigen::Matrix3Xf projectionDirections_matrix;

        const Lattice_t& l = lattice;
        Eigen::Matrix3f basis;
        basis << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;
        Lattice lattice_class(basis);

        simpleMonochromaticDiffractionPatternPrediction->predictPattern(predictedPeaks_m_matrix, millerIndices_matrix, projectionDirections_matrix,
                                                                        lattice_class);

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

} // namespace xgandalf