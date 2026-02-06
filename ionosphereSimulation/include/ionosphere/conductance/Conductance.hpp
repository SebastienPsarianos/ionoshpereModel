#pragma once
#include "ionosphere/utils/Grid.hpp"

class Conductance {
  public:
    Conductance(Grid<Sigma>& sigma, Grid<DSigma>& dConductance, size_t nTh,
                size_t nPh, double sig0, Grid<ThPh>& coords, int year, int day,
                int month, int hour);
    void calculateCoefficients();

  private:
    void _calcSigma();
    void _calcSigmaDer();
    int _auroralConductance();

    void _euvConductance(double f107);

    size_t _nTh;
    size_t _nPh;

    double _utTime;
    double _ttTime;
    double _sig0;

    Grid<ThPh>& _coords;
    Grid<Sigma>& _sigma;
    Grid<HppSigma> _hppSigma;
    Grid<DSigma> _dSigma;
};
