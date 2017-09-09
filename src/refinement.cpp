#include "refinement.h"


using namespace Eigen;
using namespace std;


// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
// return: current defect
float getGradient_reciprocalPeakMatch_meanDist(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
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

            gradient(row, col) = (numerator * denominator_inv).sum();
        }
    }

    return d_norms.array().mean();
}

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
// return: current defect
float getGradient_reciprocalPeakMatch_meanSquaredDist(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
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

	return d_norms.array().mean();

	return 0;
}

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
// return: current defect
float getGradient_detectorAngleMatch(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    Matrix2Xf detectorPeakDirections = N.bottomRows(2).colwise().normalized();

    Matrix3Xf predictedPoints = B * M;
    RowVectorXf projectionDist = (predictedPoints.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
    Matrix2Xf predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionDist.array();
    RowVectorXf defect = (predictedPoints.bottomRows(2) - predictedPointsProjected).colwise().norm();
    float summedDefect = defect.sum();

    // now do the numeric differentiatiation
    float differentiatiationShift = 0.05; // 0.05A should be small enough, but not too small
    float differentiatiationShift_inv = 1 / differentiatiationShift;
    Matrix3Xf predictedPoints_shifted;
    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            predictedPoints_shifted = predictedPoints;
            predictedPoints_shifted.row(row) += differentiatiationShift * predictedPoints.row(col);

            projectionDist = (predictedPoints_shifted.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
            predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionDist.array();
            defect = (predictedPoints_shifted.bottomRows(2) - predictedPointsProjected).colwise().norm();
            float summedDefect_shifted = defect.sum();

            gradient(row, col) = (summedDefect_shifted - summedDefect) * differentiatiationShift_inv;
        }
    }

    // return summedDefect;
    return 1;
}

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void refineReciprocalBasis_meanSquaredDist(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
{
    B = ((M * M.transpose()).colPivHouseholderQr().solve(M * N.transpose())).transpose();
}