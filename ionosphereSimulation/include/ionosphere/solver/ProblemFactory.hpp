#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/Problem.hpp"

namespace Ionosphere {
Teuchos::RCP<Problem>
problemFactory(const Teuchos::ParameterList& equationParams,
               Teuchos::RCP<Coordinates> coords,
               Ionosphere::VectorRCP radCurrent,
               Ionosphere::MultiVectorRCP conductance, Ionosphere::MapRCP map);

};
