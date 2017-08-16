/*
 * Lattice.h
 *
 *  Created on: 08.04.2017
 *      Author: Yaro
 */

#ifndef LATTICE_H_
#define LATTICE_H_

#include <Eigen/Dense>
#include <iostream>

class Lattice
{
  public:
    Lattice();
    Lattice(const Eigen::Matrix3f& basis);
    Lattice(const Eigen::Vector3f& a, const Eigen::Vector3f& b, const Eigen::Vector3f& c);

    Lattice& minimize();

    inline float det() const
    {
        return basis.determinant();
    }

    inline const Eigen::Matrix3f& getBasis() const
    {
        return basis;
    }

    inline Eigen::Vector3f getBasisVectorNorms() const
    {
        return basis.colwise().norm();
    }

    Eigen::Vector3f getBasisVectorAngles_deg() const;
    Eigen::Vector3f getBasisVectorAnglesNormalized_deg() const;

    inline Lattice getReciprocalLattice() const
    {
        return Lattice(basis.transpose().inverse().eval());
    }

    friend std::ostream& operator<<(std::ostream& os, const Lattice& lattice);

    void reorder(const Eigen::Vector3f prototypeNorms, const Eigen::Vector3f prototypeAngles_deg);

  private:
    Eigen::Matrix3f basis;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

#endif /* LATTICE_H_ */
