#include "ionosphere/solver/Problem.hpp"
#include "ionosphere/solver/BoundaryCondition.hpp"
#include "ionosphere/solver/Equation.hpp"

Problem::Problem(Teuchos::RCP<Equation> eq, Teuchos::RCP<BoundaryCondition> bc)
    : _eq(eq), _bc(bc) {}

void Problem::build() {
    // Create the problem from the equation
    matrix = _eq->assembleMatrix();
    rhs = _eq->assembleRHS();
    guess = _eq->initialGuess();

    // Apply the boundary condition
    _bc->apply(matrix, rhs);

    matrix->fillComplete();
}
