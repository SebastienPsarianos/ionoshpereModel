#include "Tpetra_Vector_decl.hpp"
#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/postProcessing/EField.hpp"
#include "ionosphere/solver/GsSolver.hpp"
#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/utils/Coordinates.hpp"
#include "ionosphere/utils/Grid.hpp"
#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Tpetra_Core.hpp>
#include <fstream>
#include <iostream>

double SIG0 = 1000;

int YEAR = 2026;
int MONTH = 2;
int DAY = 5;
int HOUR = 5;

int KP = 1;
int F107 = 120;

int main(int argc, char* argv[]) {
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();

    if (argc < 2) {
        std::cerr << "Missing filename for Radial Current data" << std::endl;
        return -1;
    }

    int nTh = -1;
    int nPh = -1;
    double dTh = -1;
    double dPh = -1;

    if (comm->getRank() == 0) {
        LegacyMHDConversion::getGridSize(argv[1], &nTh, &nPh);
        if (nTh <= 0 || nPh <= 0) {
            std::cerr << "Error reading file ending execution" << std::endl;
            return 1;
        }
        std::cout << "Sucessfully able to open radial current data, (" << nTh
                  << ", " << nPh << ") grid detected" << std::endl;
    }

    Grid<GeoSph> coords = Grid<GeoSph>(nTh, nPh);
    Grid<double> radCurrent = Grid<double>(nTh, nPh, 0.0);
    Grid<HppSigma> conductance = Grid<HppSigma>(nTh, nPh);
    Grid<DSigma> dConductance = Grid<DSigma>(nTh, nPh);
    Grid<double> potential = Grid<double>(nTh, nPh, 0.0);
    Grid<GeoSph> eField = Grid<GeoSph>(nTh, nPh);

    comm->broadcast(0, sizeof(int), (char*)&nTh);
    comm->broadcast(0, sizeof(int), (char*)&nPh);

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

    LegacyMHDConversion::processLegacyOutput(
        std::string(argv[1]), dTh, dPh, coords2, sourceTerm, comm, nTh, nPh);

    EmpConductance(conductance2, coords2, nTh, nPh, SIG0, YEAR, DAY, MONTH,
                   HOUR, map)
        .computeConductance(KP, F107);

    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        nTh * nPh, (comm->getRank() == 0 ? nTh * nPh : 0), 0, comm));
    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootSigma = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 3));

    Tpetra::Export<int, long long> exporter(map, rootMap);

    rootCoords->doExport(*coords2, exporter, Tpetra::INSERT);
    rootSigma->doExport(*conductance2, exporter, Tpetra::INSERT);

    if (comm->getRank() == 0) {
        auto sigPVals = rootSigma->getDataNonConst(0);
        auto sigHVals = rootSigma->getDataNonConst(1);
        auto sigPaVals = rootSigma->getDataNonConst(2);
        auto thVals = rootCoords->getDataNonConst(0);
        auto phVals = rootCoords->getDataNonConst(1);

        for (int i = 0; i < thVals.size(); i++) {
            int th = i % nTh;
            int ph = i / nTh;

            coords(th, ph).theta = thVals[i];
            coords(th, ph).phi = phVals[i];
            conductance(th, ph).pederson = sigPVals[i];
            conductance(th, ph).parallel = sigPaVals[i];
            conductance(th, ph).hall = sigHVals[i];
        }
        std::ofstream test("../data/test.txt");
        conductance.printWithCoords(test, coords);
        test.close();
    }

    // INPR
    TlSolver test = TlSolver(nTh, nPh, dTh, dPh, map);

    return 0;

    // TODO: Proper post-processing, including J and others
    calculateEField(eField, potential, coords);

    std::ofstream potentialOutput("../data/solvedPotential.txt");
    std::ofstream eFieldOutput("../data/solvedEField.txt");

    if (potentialOutput.is_open()) {
        potential.printWithCoords(potentialOutput, coords);
    } else {
        std::cerr << "Unable to open potential file" << std::endl;
    }

    if (eFieldOutput.is_open()) {
        eField.printWithCoords(eFieldOutput, coords);
    } else {
        std::cerr << "Unable to open eField file" << std::endl;
    }
}
