#pragma once
#include "ionosphere/utils/Grid.hpp"

enum Algorithm { GAUSS_SEIDEL, JACOBI, SOR };

class Solver {
  public:
    Solver(Grid<double>& potential, size_t nTh, size_t nPh, Grid<ThPh>& coords,
           Grid<Coeff>& kappa, Grid<double>& radCurrent,
           Grid<Sigma>& conductance, Grid<DSigma>& dConductance,
           Algorithm algorithm);
    void calculatePotential();

  private:
    size_t _nTh;
    size_t _nPh;
    double _dTh;
    double _dPh;
    double _dTh2;
    double _dPh2;

    Grid<double>& _radCurrent;
    Grid<ThPh>& _coords;
    Grid<Coeff>& _kappa;
    Grid<double>& _potential;
    Grid<Sigma>& _conductance;
    Grid<DSigma>& _dConductance;
    Algorithm _algorithm;

    double _gaussSeidelFormula(size_t th, size_t ph);
    double _calculateResidual();
    void _calcCoeff();
};
