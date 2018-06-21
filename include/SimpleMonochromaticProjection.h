#pragma once

#include "ReciprocalToRealProjection.h"
#include <Eigen/Dense>

class SimpleMonochromaticProjection : public ReciprocalToRealProjection
{
  public:
	  SimpleMonochromaticProjection(const ExperimentSettings& experimentSettings);

    void project(Eigen::Matrix2Xf& projectedPeaks, const Eigen::Matrix3Xf& reciprocalPeaks);

  private:
	  float reciprocalLambda_1A;
};
