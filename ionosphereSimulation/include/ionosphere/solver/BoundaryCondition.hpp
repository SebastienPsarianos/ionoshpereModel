#pragma once
#include "ionosphere/TrilinosAliases.hpp"

class BoundaryCondition {
  public:
    virtual ~BoundaryCondition() = default;
    virtual void apply(Ionosphere::MatrixRCP matrix,
                       Ionosphere::VectorRCP rhs) = 0;
};
