/*
 * samplePointsFiltering.h
 *
 *  Created on: 17.06.2017
 *      Author: Yaro
 */

#ifndef SAMPLEPOINTSFILTERING_H_
#define SAMPLEPOINTSFILTERING_H_

#include <Eigen/Dense>
#include "BadInputException.h"

void filterSamplePointsForNorm(Eigen::Matrix3Xf& samplePoints, const Eigen::ArrayXf& allowedNorms, float allowedTolerance);



#endif /* SAMPLEPOINTSFILTERING_H_ */
