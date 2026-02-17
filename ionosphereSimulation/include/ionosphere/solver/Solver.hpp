#pragma once
#include <cstddef>

class Solver {
  public:
    Solver(size_t nTh, size_t nPh) {
        _nTh = nTh;
        _nPh = nPh;
    }

    virtual void calculatePotential() = 0;

  protected:
    size_t _nTh;
    size_t _nPh;
};
