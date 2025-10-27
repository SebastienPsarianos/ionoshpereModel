#pragma once
#include "Grid.hpp"

class PotentialSolver {
  public:
    PotentialSolver(ThreeGrid& Sigma) : _Sigma(Sigma) {}

    int solve();

  private:
    ThreeGrid& _Sigma;
};
