#include "adaptions/crystfel/Lattice.h"

#include "Lattice.h"
#include <Eigen/Dense>

extern "C" void reorderLattice(const Lattice_t* prototype, Lattice_t* lattice)
{
    Eigen::Matrix3f prototypeBasis;
    prototypeBasis << prototype->ax, prototype->bx, prototype->cx, prototype->ay, prototype->by, prototype->cy, prototype->az, prototype->bz, prototype->cz;
    Lattice prototypeLattice(prototypeBasis);

    Eigen::Matrix3f basis;
    basis << lattice->ax, lattice->bx, lattice->cx, lattice->ay, lattice->by, lattice->cy, lattice->az, lattice->bz, lattice->cz;
    Lattice latticeWrapper(basis);

    latticeWrapper.reorder(prototypeLattice.getBasisVectorNorms(), prototypeLattice.getBasisVectorAngles_deg());
    basis = latticeWrapper.getBasis();

    lattice->ax = basis(0, 0);
    lattice->ay = basis(1, 0);
    lattice->az = basis(2, 0);
    lattice->bx = basis(0, 1);
    lattice->by = basis(1, 1);
    lattice->bz = basis(2, 1);
    lattice->cx = basis(0, 2);
    lattice->cy = basis(1, 2);
    lattice->cz = basis(2, 2);
}


void reduceLattice(Lattice_t* lattice, LatticeTransform_t* appliedReductionTransform)
{
    Eigen::Matrix3f basis, reducedBasis;
    basis << lattice->ax, lattice->bx, lattice->cx, lattice->ay, lattice->by, lattice->cy, lattice->az, lattice->bz, lattice->cz;
    Lattice latticeWrapper(basis);

    reducedBasis = latticeWrapper.minimize().getBasis();

    lattice->ax = reducedBasis(0, 0);
    lattice->ay = reducedBasis(1, 0);
    lattice->az = reducedBasis(2, 0);
    lattice->bx = reducedBasis(0, 1);
    lattice->by = reducedBasis(1, 1);
    lattice->bz = reducedBasis(2, 1);
    lattice->cx = reducedBasis(0, 2);
    lattice->cy = reducedBasis(1, 2);
    lattice->cz = reducedBasis(2, 2);

    Eigen::Matrix3f transformMatrix = basis.inverse() * reducedBasis;

    appliedReductionTransform->matrixElement_0_0 = transformMatrix(0, 0);
    appliedReductionTransform->matrixElement_1_0 = transformMatrix(1, 0);
    appliedReductionTransform->matrixElement_2_0 = transformMatrix(2, 0);
    appliedReductionTransform->matrixElement_0_1 = transformMatrix(0, 1);
    appliedReductionTransform->matrixElement_1_1 = transformMatrix(1, 1);
    appliedReductionTransform->matrixElement_2_1 = transformMatrix(2, 1);
    appliedReductionTransform->matrixElement_0_2 = transformMatrix(0, 2);
    appliedReductionTransform->matrixElement_1_2 = transformMatrix(1, 2);
    appliedReductionTransform->matrixElement_2_2 = transformMatrix(2, 2);
}

void restoreLattice(Lattice_t* lattice, LatticeTransform_t* appliedReductionTransform)
{
    Eigen::Matrix3f transformMatrix;

    transformMatrix(0, 0) = appliedReductionTransform->matrixElement_0_0;
    transformMatrix(1, 0) = appliedReductionTransform->matrixElement_1_0;
    transformMatrix(2, 0) = appliedReductionTransform->matrixElement_2_0;
    transformMatrix(0, 1) = appliedReductionTransform->matrixElement_0_1;
    transformMatrix(1, 1) = appliedReductionTransform->matrixElement_1_1;
    transformMatrix(2, 1) = appliedReductionTransform->matrixElement_2_1;
    transformMatrix(0, 2) = appliedReductionTransform->matrixElement_0_2;
    transformMatrix(1, 2) = appliedReductionTransform->matrixElement_1_2;
    transformMatrix(2, 2) = appliedReductionTransform->matrixElement_2_2;

    Eigen::Matrix3f basis, reducedBasis;
    reducedBasis << lattice->ax, lattice->bx, lattice->cx, lattice->ay, lattice->by, lattice->cy, lattice->az, lattice->bz, lattice->cz;

    basis = reducedBasis * transformMatrix.inverse();

    lattice->ax = basis(0, 0);
    lattice->ay = basis(1, 0);
    lattice->az = basis(2, 0);
    lattice->bx = basis(0, 1);
    lattice->by = basis(1, 1);
    lattice->bz = basis(2, 1);
    lattice->cx = basis(0, 2);
    lattice->cy = basis(1, 2);
    lattice->cz = basis(2, 2);
}

