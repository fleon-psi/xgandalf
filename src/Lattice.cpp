/*
 * Lattice.cpp
 *
 *  Created on: 08.04.2017
 *      Author: Yaro
 */

#include <Lattice.h>
#include <algorithm>
#include <assert.h>
#include <inttypes.h>
#include <limits>
#include <math.h>

#define M_1_PI 0.318309886183790671538

using namespace std;
using namespace Eigen;

static inline void minimize2DLattice(Matrix3f& basis);
static inline void sortTwoColumns_ascending(Matrix3f& basis, Vector3f& squaredNorms, float colNumber1, float colNumber2);
static inline void sortColumnsByNorm_ascending(Matrix3f& basis);
static inline void getMinA(Matrix3f& basis, Vector3f& minA, float& minALengthSquared);

Lattice::Lattice()
{
}

Lattice::Lattice(const Matrix3f& basis)
    : basis(basis)
{
}
Lattice::Lattice(const Vector3f& a, const Vector3f& b, const Vector3f& c)
{
    basis << a, b, c;
}

static inline void minimize2DLattice(Matrix3f& basis)
{
    if (basis.col(0).squaredNorm() < basis.col(1).squaredNorm())
    { // will always be true in the current arrangement!
        basis.col(0).swap(basis.col(1));
    }

    do
    {
        float r = round(basis.col(0).dot(basis.col(1)) / basis.col(1).squaredNorm());
        Vector3f tmp = basis.col(0) - r * basis.col(1);

        basis.col(0) = basis.col(1);
        basis.col(1) = tmp;
    } while (basis.col(1).squaredNorm() < basis.col(0).squaredNorm());
}

static inline void sortTwoColumns_ascending(Matrix3f& basis, Vector3f& squaredNorms, float colNumber1, float colNumber2)
{
    if (squaredNorms[colNumber1] > squaredNorms[colNumber2])
    {
        basis.col(colNumber1).swap(basis.col(colNumber2));

        float tmp = squaredNorms[colNumber1];
        squaredNorms[colNumber1] = squaredNorms[colNumber2];
        squaredNorms[colNumber2] = tmp;
    }
}

static inline void sortColumnsByNorm_ascending(Matrix3f& basis)
{
    // simple bubble sort implementation
    Vector3f squaredNorms = basis.colwise().squaredNorm();
    sortTwoColumns_ascending(basis, squaredNorms, 0, 1);
    sortTwoColumns_ascending(basis, squaredNorms, 1, 2);
    sortTwoColumns_ascending(basis, squaredNorms, 0, 1);
}

static inline void getMinA(Matrix3f& basis, Vector3f& minA, float& minALengthSquared)
{
    const Vector3f& b1 = basis.col(0);
    const Vector3f& b2 = basis.col(1);
    const Vector3f& b3 = basis.col(2);

    float b232 = b2.dot(b3) / b2.squaredNorm();
    float b122 = b1.dot(b2) / b2.squaredNorm();
    float b131 = b1.dot(b3) / b1.squaredNorm();
    float b121 = b1.dot(b2) / b1.squaredNorm();

    float y2 = -(b232 - b122 * b131) / (1 - b121 * b122);
    float y1 = -(b131 - b121 * b232) / (1 - b121 * b122);

    minALengthSquared = numeric_limits<float>::max();
    for (float i_1 = -1; i_1 <= 1; i_1++)
    {
        for (float i_2 = -1; i_2 <= 1; i_2++)
        {
            float x1 = round(y1 + i_1);
            float x2 = round(y2 + i_2);
            if (abs(x1 - y1) <= 1 & abs(x2 - y2) <= 1)
            {
                Vector3f a = b3 + x2 * b2 + x1 * b1;
                float aLengthSquared = a.squaredNorm();

                if (aLengthSquared < minALengthSquared)
                {
                    minA = a;
                    minALengthSquared = aLengthSquared;
                }
            }
        }
    }
}

// algorithm implemented after http://www.csie.nuk.edu.tw/~cychen/Lattices/A%203-Dimensional%20Lattice%20Reduction%20Algorithm.pdf
Lattice& Lattice::minimize()
{
    assert(abs(basis.determinant()) >= 1e-10); // nonsingular (value for normal crystal... can bee too high!)

    Vector3f minA(0, 0, 0); // initialization not required, but if not done, compiler issues warning

    bool terminationFlag = false;
    while (terminationFlag == false)
    {
        sortColumnsByNorm_ascending(basis);

        minimize2DLattice(basis);

        float minALengthSquared;
        getMinA(basis, minA, minALengthSquared);

        if (minALengthSquared >= basis.col(2).squaredNorm())
        {
            terminationFlag = true;
        }
        else
        {
            basis.col(2) = minA;
        }
    }

    sortColumnsByNorm_ascending(basis);

    return *this;
}

Vector3f Lattice::getBasisVectorAngles() const
{
    const Vector3f& a = basis.col(0);
    const Vector3f& b = basis.col(1);
    const Vector3f& c = basis.col(2);

    Vector3f angles;
    angles[0] = atan2(b.cross(c).norm(), b.dot(c)) * (M_1_PI * 180);
    angles[1] = atan2(a.cross(c).norm(), a.dot(c)) * (M_1_PI * 180);
    angles[2] = atan2(a.cross(b).norm(), a.dot(b)) * (M_1_PI * 180);

    return angles;
}

Vector3f Lattice::getBasisVectorAnglesNormalized() const
{
    const Vector3f& a = basis.col(0);
    const Vector3f& b = basis.col(1);
    const Vector3f& c = basis.col(2);

    Vector3f angles;
    angles[0] = atan2(b.cross(c).norm(), abs(b.dot(c))) * (M_1_PI * 180);
    angles[1] = atan2(a.cross(c).norm(), abs(a.dot(c))) * (M_1_PI * 180);
    angles[2] = atan2(a.cross(b).norm(), abs(a.dot(b))) * (M_1_PI * 180);

    return angles;
}

std::ostream& operator<<(std::ostream& os, const Lattice& lattice)
{
    os << lattice.basis;
    return os;
}
