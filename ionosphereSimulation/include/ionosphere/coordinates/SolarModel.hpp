#pragma once
#include "ionosphere/coordinates/Grid.hpp"

#define SUNFREQ .0172019715

class SolarModel {
  public:
    SolarModel(int year, int month, int day, int hour);
    double computeZenith(GeoGeo coords) const;
    GeoGeo computeSubSolar() const;

  private:
    void _computeSunPosition(double& alpha, double& declination) const;
    double _utTime;
    double _ttTime;
};
