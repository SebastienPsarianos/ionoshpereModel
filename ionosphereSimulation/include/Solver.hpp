#pragma once

#include "Grid.hpp"
#include <memory>

enum Algorithm { GAUSS_SEIDEL, JACOBI, SOR };

class Solver {
  public:
    Solver(size_t nTh, size_t nPh, std::shared_ptr<GridSet<Coeff>> kappa,
           std::shared_ptr<GridSet<Ang>> coords,
           std::shared_ptr<Grid> radCurrent, Algorithm algorithm);
    std::shared_ptr<Grid> calculatePotential();

  private:
    size_t nTh;
    size_t nPh;
    double dTh;
    double dPh;
    double dTh2;
    double dPh2;
    std::shared_ptr<Grid> radCurrent;
    std::shared_ptr<GridSet<Ang>> coords;
    std::shared_ptr<GridSet<Coeff>> kappa;
    std::shared_ptr<Grid> potential;
    Algorithm algorithm;
    Grid previousIteration;

    double _gaussSeidelFormula(size_t th, size_t ph);
    double _calculateResidual();
};
