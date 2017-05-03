/*
 * ReciprocalSpaceTranslation.cpp
 *
 *  Created on: 01.05.2017
 *      Author: Yaro
 */

#include <DetectorToReciprocalSpaceTransform.h>

using namespace Eigen;
using namespace std;

DetectorToReciprocalSpaceTransform::DetectorToReciprocalSpaceTransform(const ExperimentSettings& experimentSettings)
{
    reciprocal_lambda_1A = experimentSettings.getReciprocalLambda_1A();
    detectorDistance_m = experimentSettings.getDetectorDistance_m();
}

void DetectorToReciprocalSpaceTransform::computeReciprocalPeaksFromDetectorPeaks(Matrix3Xf& reciprocalPeaks_A, const Matrix2Xf& detectorPeaks_m)
{
    Matrix3Xf backprojectionDirectionVectors(3, detectorPeaks_m.cols());
    backprojectionDirectionVectors.row(0).setConstant(detectorDistance_m);
    backprojectionDirectionVectors.row(1) = -1*detectorPeaks_m.row(0);  // detector x-coordinate is -y coordinate in reciprocal space
    backprojectionDirectionVectors.row(2) = detectorPeaks_m.row(1);

    reciprocalPeaks_A = backprojectionDirectionVectors.colwise().normalized() * reciprocal_lambda_1A;
    reciprocalPeaks_A.row(0) -= RowVectorXf::Constant(reciprocalPeaks_A.cols(), reciprocal_lambda_1A);
}
