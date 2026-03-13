#pragma once
#include "ionosphere/coordinates/CoordinateTypes.hpp"

#include <Eigen/Dense>

class DipoleModel {
  public:
    DipoleModel();
    Ionosphere::MagSph geoCentricToDipole(Ionosphere::GeoSph geoCoords) const;
    Ionosphere::GeoSph dipoleToGeoCentric(Ionosphere::MagSph magCoords) const;

  private:
    Eigen::Matrix3d _rotationMatrix;
};
