#pragma once
#include "ionosphere/coordinates/Grid.hpp"

#include <Eigen/Dense>

class DipoleModel {
  public:
    DipoleModel();
    MagSph geoCentricToDipole(GeoSph geoCoords) const;
    GeoSph dipoleToGeoCentric(MagSph magCoords) const;

  private:
    Eigen::Matrix3d _rotationMatrix;
};
