#pragma once
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/solver/BoundaryCondition.hpp"
#include "ionosphere/solver/Equation.hpp"

class Problem {
  public:
    Problem(Teuchos::RCP<Equation> eq, Teuchos::RCP<BoundaryCondition> bc);

    void build();

    Ionosphere::MatrixRCP matrix;
    Ionosphere::VectorRCP rhs;
    Ionosphere::VectorRCP guess;

  private:
    Teuchos::RCP<Equation> _eq;
    Teuchos::RCP<BoundaryCondition> _bc;
};
