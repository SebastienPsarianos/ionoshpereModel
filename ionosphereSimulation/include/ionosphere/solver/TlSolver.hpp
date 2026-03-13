#pragma once

#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/Solver.hpp"

class TlSolver : Solver {
  public:
    TlSolver(Teuchos::RCP<Coordinates> coords,
             Ionosphere::MultiVectorRCP conductance,
             Ionosphere::VectorRCP sourceTerm, Ionosphere::MapRCP map);

    Ionosphere::VectorRCP calculatePotential();

  private:
    Ionosphere::MultiVectorRCP _calculateCoefficients();
    // TODO: Might not need this? Ionosphere::VectorRCP _buildSourceVector();
    Ionosphere::MatrixRCP _buildGrid(Ionosphere::MultiVectorRCP coefficients);
    Ionosphere::MultiVectorRCP _gatherPoleData();

    Teuchos::RCP<Coordinates> _coords;
    Ionosphere::MultiVectorRCP _conductance;
    Ionosphere::VectorRCP _sourceTerm;
    Ionosphere::MapRCP _map;
};
