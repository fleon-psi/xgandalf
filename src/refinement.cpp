#include "refinement.h"


using namespace Eigen;
using namespace std;


// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void getGradient_reciprocalPeakMatch_meanDist(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    Matrix3Xf d = B * M - N;
    RowVectorXf d_norms = d.colwise().norm();

    Array<float, 1, Dynamic> denominator_inv = 1 / (d_norms.array());

    for (int i = 0; i < denominator_inv.size(); ++i)
    {
        if (!isfinite(denominator_inv[i]))
        {
            denominator_inv = 1e-30;
        }
    }

    Array<float, 1, Dynamic> numerator;
    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            numerator = d.row(row).array() * M.row(col).array();

            gradient(row, col) = (numerator * denominator_inv).mean();
        }
    }
}


// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
// TODO: Function can be made faster easily... nevertheless double is really necessary...
void getGradient_detectorAngleMatch(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    Matrix2Xd detectorPeakDirections = N.bottomRows(2).colwise().normalized().cast<double>();

    Matrix3Xd Md = M.cast<double>();

    Matrix3Xd predictedPoints = B.cast<double>() * Md.cast<double>();
    RowVectorXd projectionNorms = (predictedPoints.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
    Matrix2Xd predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionNorms.array();
    RowVectorXd defect = (predictedPoints.bottomRows(2) - predictedPointsProjected).colwise().norm();
    double meanDefect = defect.mean();

    // now do the numeric differentiatiation
    double differentiatiationShift = 1e-8; // should be small enough, but not too small
    double differentiatiationShift_inv = 1 / differentiatiationShift;
    Matrix3Xd predictedPoints_shifted;
    for (int row = 1; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            predictedPoints_shifted = predictedPoints;
            predictedPoints_shifted.row(row) += differentiatiationShift * Md.row(col);


            projectionNorms = (predictedPoints_shifted.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
            predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionNorms.array();
            defect = (predictedPoints_shifted.bottomRows(2) - predictedPointsProjected).colwise().norm();
            double meanDefect_shifted = defect.mean();

            gradient(row, col) = (meanDefect_shifted - meanDefect) * differentiatiationShift_inv;
        }
    }

    // shifts in x direction do not change the angle
    gradient(0, 0) = 0;
    gradient(0, 1) = 0;
    gradient(0, 2) = 0;
}

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void getGradient_reciprocalPeakMatch_meanSquaredDist(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    Matrix3Xf d = B * M - N;
    RowVectorXf d_norms = d.colwise().norm();

    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            gradient(row, col) = (d.row(row).array() * M.row(col).array()).sum();
        }
    }
}

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void refineReciprocalBasis_meanSquaredDist(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    B = M.transpose().colPivHouseholderQr().solve(N.transpose()).transpose();
}

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void refineReciprocalBasis_meanDist_peaksAndAngle(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    Matrix3f gradient;
    Matrix3f summedGradient;

    float reciprocalPeakDistWeight = 1;
    float reciprocalPeakAngleWeight = 2;

    float stepLength = B.maxCoeff() * 0.003;
    for (int i = 0; i < 150; i++)
    {
        getGradient_reciprocalPeakMatch_meanDist(gradient, B, M, N);
        summedGradient = gradient * reciprocalPeakDistWeight;
        getGradient_detectorAngleMatch(gradient, B, M, N);
        summedGradient += gradient * reciprocalPeakAngleWeight;
        float maxCoeff = summedGradient.cwiseAbs().maxCoeff();
        if (maxCoeff < 1e-10)
        {
            break;
        }
        summedGradient /= maxCoeff;
        B = B - stepLength * summedGradient;

        if (i >= 75)
        {
            stepLength *= 0.93;
        }
    }
}