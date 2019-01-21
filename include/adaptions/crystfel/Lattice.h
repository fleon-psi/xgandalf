/* 
 * Lattice.h
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