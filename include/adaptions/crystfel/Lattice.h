#ifndef ADAPTIONS_CRYSTFEL_LATTICE_H
#define ADAPTIONS_CRYSTFEL_LATTICE_H

typedef struct
{
    float ax;
    float ay;
    float az;

    float bx;
    float by;
    float bz;

    float cx;
    float cy;
    float cz;
} Lattice_t;

typedef struct
{
    float matrixElement_0_0;
    float matrixElement_1_0;
    float matrixElement_2_0;
    float matrixElement_0_1;
    float matrixElement_1_1;
    float matrixElement_2_1;
    float matrixElement_0_2;
    float matrixElement_1_2;
    float matrixElement_2_2;
} LatticeTransform_t;

#ifdef __cplusplus
extern "C" {
#endif

void reorderLattice(const Lattice_t* prototype, Lattice_t* lattice);
void reduceLattice(Lattice_t* lattice, LatticeTransform_t* appliedReductionTransform);
void restoreLattice(Lattice_t* lattice, LatticeTransform_t* appliedReductionTransform);

#ifdef __cplusplus
}
#endif

#endif