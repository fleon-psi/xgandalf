#include "refinement.h"


using namespace Eigen;
using namespace std;


// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void getGradient(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
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
}