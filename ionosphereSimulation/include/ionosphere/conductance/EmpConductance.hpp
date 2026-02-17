#pragma once
#include "ionosphere/utils/Grid.hpp"

class EmpConductance {
  public:
    EmpConductance(size_t nTh, size_t nPh, double sig0, int year, int day,
                   int month, int hour);
    int computeConductance(Grid<Sigma>& sigma, Grid<DSigma>& dConductance,
                           Grid<GeoSph>& coords);

  private:
    void _calcSigma(Grid<Sigma>& sigma, Grid<GeoSph>& coords);
    void _calcSigmaDer(Grid<DSigma>& dConductance, Grid<Sigma>& sigma,
                       Grid<GeoSph>& coords);

    int _computeAuroralConductance(int kp, Grid<GeoSph>& coords);
    void _computeEuvConductance(double f107, Grid<GeoSph>& coords);
    void _computeHppConductance();

    size_t _nTh;
    size_t _nPh;

    double _utTime;
    double _ttTime;
    double _sig0;

    Grid<HppSigma> _hppSigma;
    Grid<HppSigma> _euvConductance;
    Grid<HppSigma> _auroralConductance;
};
