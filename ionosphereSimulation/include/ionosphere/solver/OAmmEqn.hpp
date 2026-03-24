#pragma once

#include "Teuchos_ParameterList.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/Equation.hpp"

class OAmmEqn : public Equation {
  public:
    OAmmEqn(const Teuchos::ParameterList& equationParams,
            Teuchos::RCP<Coordinates> coords,
            Ionosphere::MultiVectorRCP conductance,
            Ionosphere::VectorRCP radCurrent, Ionosphere::MapRCP map);

    Ionosphere::VectorRCP assembleRHS() override;
    Ionosphere::MatrixRCP assembleMatrix() override;
    Ionosphere::VectorRCP initialGuess() override;

  private:
    Ionosphere::MultiVectorRCP _gatherPoleData();
    Ionosphere::MultiVectorRCP _calculateCoefficients();

    Teuchos::RCP<Coordinates> _coords;
    Ionosphere::MultiVectorRCP _conductance;
    Ionosphere::MultiVectorRCP _coefficients;
    Ionosphere::VectorRCP _radCurrent;
    Ionosphere::MapRCP _map;

    const Ionosphere::Scalar _radiusIonosphere2;
};
