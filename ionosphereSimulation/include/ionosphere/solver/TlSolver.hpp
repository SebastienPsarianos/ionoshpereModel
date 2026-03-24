#pragma once

#include "Teuchos_RCPDecl.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/solver/Problem.hpp"

namespace Ionosphere {
Ionosphere::VectorRCP calculatePotential(Teuchos::RCP<Problem> problem);
}
