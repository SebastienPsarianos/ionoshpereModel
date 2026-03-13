#pragma once

#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/Solver.hpp"

class TlSolver : Solver {
  public:
    TlSolver(Teuchos::RCP<Coordinates> coords,
             Ionosphere::MultiVectorRCP conductance,
             Ionosphere::VectorRCP sourceTerm, Ionosphere::MapRCP map);

    Ionosphere::VectorRCP calculatePotential();

    /// Build the stiffness matrix on the current map (exposed for repartitioning)
    Ionosphere::MatrixRCP buildMatrix();

  private:
    Ionosphere::MultiVectorRCP _calculateCoefficients();
    Ionosphere::MatrixRCP _buildGrid(Ionosphere::MultiVectorRCP coefficients);
    Ionosphere::MultiVectorRCP _gatherPoleData();

    Teuchos::RCP<Coordinates> _coords;
    Ionosphere::MultiVectorRCP _conductance;
    Ionosphere::VectorRCP _sourceTerm;
    Ionosphere::MapRCP _map;
};
