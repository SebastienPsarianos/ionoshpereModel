#pragma once
#include "Grid.hpp"

enum Algorithm { GAUSS_SEIDEL, JACOBI, SOR };

class Solver {
  public:
    Solver(Grid<double>& potential, size_t nTh, size_t nPh, Grid<Coeff>& kappa,
           Grid<ThPh>& coords, Grid<double>& radCurrent, Algorithm algorithm);
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
    Algorithm _algorithm;

    double _gaussSeidelFormula(size_t th, size_t ph);
    double _calculateResidual();
};
