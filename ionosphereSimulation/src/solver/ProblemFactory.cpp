#include "ionosphere/solver/ProblemFactory.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "ionosphere/solver/BoundaryCondition.hpp"
#include "ionosphere/solver/Equation.hpp"
#include "ionosphere/solver/GhostRingBC.hpp"
#include "ionosphere/solver/OAmmEqn.hpp"
#include "ionosphere/solver/PoleDivThmBC.hpp"
#include "ionosphere/solver/Problem.hpp"
#include "ionosphere/solver/SolverTypes.hpp"

Teuchos::RCP<Problem>
Ionosphere::problemFactory(const Teuchos::ParameterList& equationParams,
                           Teuchos::RCP<Coordinates> coords,
                           VectorRCP radCurrent, MultiVectorRCP conductance,
                           MapRCP map) {

    Teuchos::RCP<Equation> eq;
    Teuchos::RCP<BoundaryCondition> bc;

    EqnOption equation =
        equationTypeFromString(equationParams.get<std::string>("equation"));
    BCOption boundaryCondition =
        bcTypeFromString(equationParams.get<std::string>("boundarycondition"));

    switch (equation) {
    // Here different potential models can be implemented as an
    // Equation instance
    case OAMM:
        eq = Teuchos::rcp(
            new OAmmEqn(equationParams, coords, conductance, radCurrent, map));
        break;
    }

    switch (boundaryCondition) {
    // Here different boundary conditions can be
    // implemented as instances of BoundaryCondition
    case GHOSTRING:
        bc = Teuchos::rcp<BoundaryCondition>(new GhostRingBC());
        break;
    case POLEDIVTHM:
        bc = Teuchos::rcp<BoundaryCondition>(
            new PoleDivThmBC(equationParams, coords, conductance, map));
        break;
    }

    return Teuchos::rcp<Problem>(new Problem(eq, bc));
}
