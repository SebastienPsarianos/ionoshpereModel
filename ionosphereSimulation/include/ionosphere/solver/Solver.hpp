#pragma once
#include "ionosphere/TrilinosAliases.hpp"

class Solver {
    void calculatePotential(Ionosphere::MultiVectorRCP conductance,
                            Ionosphere::MultiVectorRCP coords,
                            Ionosphere::VectorRCP sourceTerm);
};
