#ifndef REFINEMENT_H_
#define REFINEMENT_H_

#include <Eigen/Dense>

// B: reciprocal basis
// M: miller indices (reciprocal)
// N: reciprocal peaks
void getGradient(Eigen::Matrix3f& gradient, const Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N);




#endif /* REFINEMENT_H_ */