#pragma once
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/CoordinateTypes.hpp"

#include <Eigen/Dense>

class DipoleModel {
  public:
    DipoleModel(Ionosphere::CommRCP comm, int year, int month, int day,
                double hour);
    Ionosphere::MagSph geoCentricToDipole(Ionosphere::GeoSph geoCoords) const;
    Ionosphere::GeoSph dipoleToGeoCentric(Ionosphere::MagSph magCoords) const;

  private:
    Eigen::Matrix3d _rotationMatrix;
};
