#include "SimpleProjection.h"

using namespace std;
using namespace Eigen;

SimpleProjection::SimpleProjection(const ExperimentSettings& experimentSettings)
    : ReciprocalToRealProjection(experimentSettings)
{
    reciprocalLambda_1A = experimentSettings.getReciprocalLambda_1A();
}


void SimpleProjection::project(Matrix2Xf& projectedPeaks, const Matrix3Xf& reciprocalPeaks)
{
    RowVectorXf yzSquaredNorms = reciprocalPeaks.bottomRows(2).colwise().squaredNorm();
    // RowVectorXf rayOriginsX = (reciprocalPeaks.row(0) + (yzSquaredNorms.array() / reciprocalPeaks.row(0).array()).matrix()) / 2;
    RowVectorXf rayOriginsX;
    rayOriginsX.setConstant(yzSquaredNorms.size(), -1 * reciprocalLambda_1A);

    projectedPeaks =
        reciprocalPeaks.bottomRows(2).array().rowwise() / (reciprocalPeaks.row(0) - rayOriginsX).array() * experimentSettings.getDetectorDistance_m();
}