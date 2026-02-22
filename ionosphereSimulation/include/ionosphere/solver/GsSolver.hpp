#pragma once
#include "ionosphere/solver/Solver.hpp"
#include "ionosphere/utils/Grid.hpp"

class GsSolver {
  public:
    GsSolver(Grid<double>& potential, size_t nTh, size_t nPh,
             Grid<GeoSph>& coords, Grid<double>& radCurrent,
             Grid<Sigma>& conductance, Grid<DSigma>& dConductance);

    void calculatePotential();

  private:
    double _dTh;
    double _dPh;
    double _dTh2;
    double _dPh2;
    double _nTh;
    double _nPh;

    Grid<double>& _radCurrent;
    Grid<GeoSph>& _coords;
    Grid<double>& _potential;
    Grid<Sigma>& _conductance;
    Grid<DSigma>& _dConductance;
    Grid<Coeff> _kappa;

    double _gaussSeidelFormula(size_t th, size_t ph);
    double _calculateResidual();
    void _calcCoeff();
};
