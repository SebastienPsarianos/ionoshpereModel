#pragma once
#include "Grid.hpp"

class Conductance {
  public:
    static void calculateConductance(GridSet<Coords>& coords,
                                     GridSet<Sigma>& sigma,
                                     GridSet<HppSigma>& hppSigma,
                                     GridSet<DSigma>& dSigma);

  private:
    static void _calcSigma(GridSet<Coords>& coords, GridSet<Sigma>& sigma,
                           GridSet<HppSigma>& hppSigma);

    static void _calcSigmaDer(GridSet<Coords>& coords, GridSet<Sigma>& sigma,
                              GridSet<DSigma>& dSigma);
};
