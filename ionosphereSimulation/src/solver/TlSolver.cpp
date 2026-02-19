#include "ionosphere/solver/TlSolver.hpp"
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

TlSolver::TlSolver(size_t nTh, size_t nPh, Tpetra::ScopeGuard scope,
                   Teuchos::RCP<const Teuchos::Comm<int>> comm)
    : Solver(nTh, nPh) {
    auto map = Teuchos::rcp(new Tpetra::Map<>(nTh, nPh, 0, comm));

    // TODO: Figure out how many non zero per row
    auto A = Teuchos::rcp(new Tpetra::CrsMatrix<>(map, 5));

    auto myElements = map->getMyGlobalIndices();

    for (size_t i = 0; i < myElements.size(); i++) {
        auto currentPoint = myElements[i];

        long long theta = currentPoint % _nTh;
        long long phi = currentPoint / _nTh;

        Teuchos::Array<long long> matrixIdcs;
        Teuchos::Array<double> vals;

        if (theta == 0) {
        }

        else if (phi == 0) {
        }

        else if (theta == _nTh - 1) {
        }

        else if (phi == _nPh - 1) {
        }

        else {
            matrixIdcs.push_back(currentPoint);
            vals.push_back(4.0);

            matrixIdcs.push_back(currentPoint - 1);
            vals.push_back(-1.0);

            matrixIdcs.push_back(currentPoint + 1);
            vals.push_back(-1.0);

            // Down
            matrixIdcs.push_back(currentPoint - _nTh);
            vals.push_back(-1.0);

            // Up
            matrixIdcs.push_back(currentPoint + _nTh);
            vals.push_back(-1.0);
        }

        Teuchos::Array<long long> indices;
        Teuchos::Array<double> values;
    }
}

void TlSolver::calculatePotential() {}
