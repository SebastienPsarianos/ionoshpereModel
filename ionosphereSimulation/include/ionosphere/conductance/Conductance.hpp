#pragma once
#include "ionosphere/utils/Grid.hpp"

class Conductance {
  public:
    Conductance(Grid<Sigma>& sigma, Grid<DSigma>& dConductance, size_t nTh,
                size_t nPh, double sig0, double sigP, double sigH,
                Grid<ThPh>& coords);
    void calculateCoefficients();

  private:
    void _calcSigma();
    void _calcSigmaDer();

    size_t _nTh;
    size_t _nPh;

    Grid<ThPh>& _coords;
    Grid<Sigma>& _sigma;
    Grid<HppSigma> _hppSigma;
    Grid<DSigma> _dSigma;
};
