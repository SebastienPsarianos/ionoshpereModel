#pragma once
#include "Coords.hpp"
#include "Grid.hpp"

class Conductance {
  public:
    void calculateConductance(std::shared_ptr<GridSet<Coords>> coords,
                              std::shared_ptr<GridSet<Sigma>> sigma,
                              std::shared_ptr<GridSet<HppSigma>> hppSigma,
                              std::shared_ptr<GridSet<DSigma>> dSigma, int nTh,
                              int nPh);

  private:
    static void _calcSigma(std::shared_ptr<GridSet<Coords>> coords,
                           std::shared_ptr<GridSet<Sigma>> sigma,
                           std::shared_ptr<GridSet<HppSigma>> hppSigma, int nTh,
                           int nPh);
    static void _calcSigmaDer(std::shared_ptr<GridSet<Coords>> coords,
                              std::shared_ptr<GridSet<Sigma>> sigma,
                              std::shared_ptr<GridSet<DSigma>> dSigma, int nTh,
                              int nPh);
};
