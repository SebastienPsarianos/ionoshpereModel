#pragma once

#define SUNFREQ .0172019715

#include "ionosphere/coordinates/CoordinateTypes.hpp"

class SolarModel {
  public:
    SolarModel(int year, int month, int day, double hour);
    double computeZenith(Ionosphere::GeoGeo coords) const;
    Ionosphere::GeoGeo computeSubSolar() const;

  private:
    void _computeSunPosition(double& alpha, double& declination) const;
    double _utTime;
    double _ttTime;
};
