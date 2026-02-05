#pragma once
#include "Grid.hpp"

class Conductance {
  public:
    Conductance(Grid<Coeff>& kappa, Grid<Sigma>& sigma, size_t nTh, size_t nPh,
                double sig0, double sigP, double sigH, Grid<ThPh>& coords);
    void calculateCoefficients();

  private:
    void _calcSigma();
    void _calcSigmaDer();
    void _calcCoeff();

    size_t _nTh;
    size_t _nPh;

    Grid<ThPh>& _coords;
    Grid<Sigma>& _sigma;
    Grid<Coeff>& _kappa;
    Grid<HppSigma> _hppSigma;
    Grid<DSigma> _dSigma;
};
