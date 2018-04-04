#pragma once

#include "ExperimentSettings.h"
#include <Eigen/Dense>


class ReciprocalToRealProjection
{
  public:
    ReciprocalToRealProjection(const ExperimentSettings& experimentSettings);
    virtual ~ReciprocalToRealProjection() = default;

    virtual void project(Eigen::Matrix2Xf& projectedPoints, const Eigen::Matrix3Xf& reciprocalPoints) = 0;

  protected:
    ExperimentSettings experimentSettings;
};
