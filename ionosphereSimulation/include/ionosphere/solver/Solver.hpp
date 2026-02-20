#pragma once
#include "ionosphere/utils/Grid.hpp"
#include <cstddef>

class Solver {
  public:
    Solver(size_t nTh, size_t nPh) {
        _nTh = nTh;
        _nPh = nPh;
    }

    virtual void calculatePotential(MultiVectorRcp radCurrent,
                                    MultiVectorRcp conductance,
                                    MultiVectorRcp coords) = 0;

  protected:
    size_t _nTh;
    size_t _nPh;
};
