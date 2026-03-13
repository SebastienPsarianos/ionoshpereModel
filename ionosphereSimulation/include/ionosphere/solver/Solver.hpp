#pragma once
#include "ionosphere/IonosphereTypes.hpp"

class Solver {
    void calculatePotential(Ionosphere::MultiVectorRCP conductance,
                            Ionosphere::MultiVectorRCP coords,
                            Ionosphere::VectorRCP sourceTerm);
};
