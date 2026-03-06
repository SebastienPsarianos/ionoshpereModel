#pragma once
#include "ionosphere/utils/Grid.hpp"

#include <Eigen/Dense>

class DipoleModel {
  public:
    DipoleModel();
    MagSph geoCentricToDipole(GeoSph geoCoords) const;
    double computeMLT(MagGeo subsolar, MagGeo observer) const;
    const Eigen::Matrix3d& rotationMatrix() const;

  private:
    Eigen::Matrix3d _rotationMatrix;
};
