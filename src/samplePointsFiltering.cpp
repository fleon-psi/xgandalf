/*
 * samplePointsFiltering.cpp
 *
 *  Created on: 17.06.2017
 *      Author: Yaro
 */

#include "samplePointsFiltering.h"

using namespace Eigen;
using namespace std;

void filterSamplePointsForNorm(Matrix3Xf& samplePoints, const ArrayXf& allowedNorms, float allowedTolerance)
{
    if (allowedNorms.size() > 3) {
        std::stringstream errStream;
        errStream << "Maximum supported number of allowedNorms is 3" << endl;
        throw BadInputException(errStream.str());
    }

    int64_t keptSamplePontsCount = 0;
    if (allowedNorms.size() == 1) {
        float allowedBordersMinimum = (allowedNorms * (1 - allowedTolerance)).square()[0];
        float allowedBordersMaximum = (allowedNorms * (1 + allowedTolerance)).square()[0];

        for (int64_t i = 0; i < samplePoints.cols(); ++i) {
            const Vector3f& samplePoint = samplePoints.col(i);
            const float squaredNorm = samplePoint.squaredNorm();
            if ((squaredNorm >= allowedBordersMinimum) & (squaredNorm <= allowedBordersMaximum)) {
                samplePoints.col(keptSamplePontsCount++) = samplePoint;
            }
        }
    } else if (allowedNorms.size() == 2) {
        Array2f allowedBordersMinima = (allowedNorms * (1 - allowedTolerance)).square();
        Array2f allowedBordersMaxima = (allowedNorms * (1 + allowedTolerance)).square();

        for (int64_t i = 0; i < samplePoints.cols(); ++i) {
            const Vector3f& samplePoint = samplePoints.col(i);
            const float squaredNorm = samplePoint.squaredNorm();
            if (((squaredNorm >= allowedBordersMinima[1]) & (squaredNorm <= allowedBordersMaxima[1])) ||
                    ((squaredNorm >= allowedBordersMinima[0]) & (squaredNorm <= allowedBordersMaxima[0]))) {
                samplePoints.col(keptSamplePontsCount++) = samplePoint;
            }
        }
    } else {
        Array3f allowedBordersMinima = (allowedNorms * (1 - allowedTolerance)).square();
        Array3f allowedBordersMaxima = (allowedNorms * (1 + allowedTolerance)).square();

        for (int64_t i = 0; i < samplePoints.cols(); ++i) {
            const Vector3f& samplePoint = samplePoints.col(i);
            const float squaredNorm = samplePoint.squaredNorm();
            if (((squaredNorm >= allowedBordersMinima[2]) & (squaredNorm <= allowedBordersMaxima[2])) ||
                    ((squaredNorm >= allowedBordersMinima[1]) & (squaredNorm <= allowedBordersMaxima[1])) ||
                    ((squaredNorm >= allowedBordersMinima[0]) & (squaredNorm <= allowedBordersMaxima[0]))) {
                samplePoints.col(keptSamplePontsCount++) = samplePoint;
            }
        }
    }

    samplePoints.conservativeResize(3, keptSamplePontsCount);

}

