#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/BoundaryCondition.hpp"

class PoleDivThmBC : public BoundaryCondition {
  public:
    PoleDivThmBC(const Teuchos::ParameterList& equationParams,
                 Teuchos::RCP<Coordinates> coords,
                 Ionosphere::MultiVectorRCP conductance,
                 Ionosphere::MapRCP map);

    void apply(Ionosphere::MatrixRCP matrix,
               Ionosphere::VectorRCP rhs) override;

  private:
    Ionosphere::MultiVectorRCP _gatherPoleData();
    void _applyToMatrix(Ionosphere::MatrixRCP matrix);
    void _applyToRHS(Ionosphere::VectorRCP rhs);

    Ionosphere::MultiVectorRCP _conductance;
    Ionosphere::MapRCP _map;
    Teuchos::RCP<Coordinates> _coords;

    const Ionosphere::Scalar _radiusIonosphere2;
};
