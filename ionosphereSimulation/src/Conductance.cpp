#include "Conductance.hpp"

Conductance::Conductance(size_t nTh, size_t nPh)
    : _conductances(3, nTh, nPh), _conductanceDerivatives(4, nTh, nPh) {}
