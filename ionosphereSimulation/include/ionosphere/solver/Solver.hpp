#pragma once
#include "ionosphere/utils/Grid.hpp"
#include <cstddef>

class Solver {
  public:
    Solver(size_t nTh, size_t nPh) {
        _nTh = nTh;
        _nPh = nPh;
    }

    void calculatePotential(MultiVectorRcp conductance, MultiVectorRcp coords,
                            VectorRcp sourceTerm);

  protected:
    size_t _nTh;
    size_t _nPh;
};
