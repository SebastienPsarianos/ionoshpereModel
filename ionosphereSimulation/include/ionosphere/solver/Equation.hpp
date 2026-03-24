#pragma once
#include "ionosphere/TrilinosAliases.hpp"

class Equation {
  public:
    virtual ~Equation() = default;
    virtual Ionosphere::VectorRCP assembleRHS() = 0;
    virtual Ionosphere::MatrixRCP assembleMatrix() = 0;
    virtual Ionosphere::VectorRCP initialGuess() = 0;
};
