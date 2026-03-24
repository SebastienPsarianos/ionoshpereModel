#pragma once

#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/solver/BoundaryCondition.hpp"

class GhostRingBC : public BoundaryCondition {
  public:
    void apply(Ionosphere::MatrixRCP matrix,
               Ionosphere::VectorRCP rhs) override;
};
