#include "ionosphere/solver/Solver.hpp"
#include <Tpetra_Core.hpp>

class TlSolver : Solver {
  public:
    TlSolver(size_t nTh, size_t nPh, Tpetra::ScopeGuard scope,
             Teuchos::RCP<const Teuchos::Comm<int>> comm);

    void calculatePotential() override;
};
