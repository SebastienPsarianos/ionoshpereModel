#pragma once
#include "Grid.hpp"

const int THTH = 0;
const int THPH = 1;
const int PHPH = 2;

const int DTHTH_TH = 0;
const int DTHPH_PH = 1;
const int DTHPH_TH = 2;
const int DPHPH_PH = 3;

class Conductance {
  public:
    Conductance(size_t nTh, size_t nPh);

    GridSet _conductances;
    GridSet _conductanceDerivatives;
};
