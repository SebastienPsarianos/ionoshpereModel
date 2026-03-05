#pragma once
#include "ionosphere/IonosphereTypes.hpp"
#include <cstddef>

class Solver {
  public:
    Solver(size_t nTh, size_t nPh) {
        _nTh = nTh;
        _nPh = nPh;
    }

    void calculatePotential(Ionosphere::MultiVectorRCP conductance,
                            Ionosphere::MultiVectorRCP coords,
                            Ionosphere::VectorRCP sourceTerm);

  protected:
    size_t _nTh;
    size_t _nPh;
};
