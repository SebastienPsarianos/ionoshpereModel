#include "Tpetra_Vector_decl.hpp"
#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/postProcessing/EField.hpp"
#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/tecplot.hpp"
#include "ionosphere/utils/Coordinates.hpp"
#include "ionosphere/utils/Grid.hpp"
#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Tpetra_Core.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

double SIG0 = 1000;
double THETA0 = 0.001;

int YEAR = 2026;
int MONTH = 2;
int DAY = 5;
int HOUR = 5;

int KP = 1;
int F107 = 120;

int main(int argc, char* argv[]) {
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    if (argc < 2) {
        std::cerr << "Missing filename for Radial Current data" << std::endl;
        return -1;
    }

    int nTh = -1;
    int nPh = -1;
    double dTh = -1;
    double dPh = -1;

    LegacyMHDConversion::getGridSize(argv[1], &nTh, &nPh, comm);
    if (nTh <= 0 || nPh <= 0) {
        throw std::runtime_error(
            "Error reading source term data, ending execution");
    }

    auto map =
        Teuchos::rcp(new Tpetra::Map<int, long long>(nTh * nPh, 0, comm));

    auto coords2 =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(map, 2));
    auto conductance2 =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(map, 3));
    auto sourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(map));
    auto potential2 =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(map));

    LegacyMHDConversion::processLegacyOutput(std::string(argv[1]), dTh, dPh,
                                             coords2, sourceTerm, comm, nTh,
                                             nPh, THETA0);

    EmpConductance(conductance2, coords2, nTh, nPh, SIG0, YEAR, DAY, MONTH,
                   HOUR, map)
        .computeConductance(KP, F107);

    // INPR
    TlSolver test = TlSolver(nTh, nPh, dTh, dPh, map);
    VectorRcp result =
        test.calculatePotential(conductance2, coords2, sourceTerm);

    exportToTecplot("../data/test3.dat", coords2, result, sourceTerm, comm, nTh,
                    nPh);
    return 0;

    // TODO: Proper post-processing, including J and others
}
