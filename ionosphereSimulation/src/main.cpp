#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/tecplot.hpp"
#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <Tpetra_Core.hpp>
#include <iostream>
#include <stdexcept>

Ionosphere::Scalar SIG0 = 1000;
Ionosphere::Scalar THETA0 = 0.05;

int YEAR = 2026;
int MONTH = 2;
int DAY = 5;
int HOUR = 5;

int KP = 6;
int F107 = 120;

int main(int argc, char* argv[]) {
    using namespace Ionosphere;
    using Teuchos::rcp;
    using Teuchos::RCP;

    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();

    /*** BEGIN PARAMETER PARSING ***/
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [-t|-m] <input_file>"
                  << std::endl;
        return -1;
    }

    bool useTecplot = true;
    std::string inputFile;

    if (argc == 2) {
        inputFile = argv[1];
    } else if (argc == 3) {
        std::string flag = argv[1];
        if (flag == "-m") {
            useTecplot = false;
        } else if (flag != "-t") {
            std::cerr << "Unknown flag: " << flag
                      << ". Use -t for TecPlot or -m for matplotlib."
                      << std::endl;
            return -1;
        }
        inputFile = argv[2];
    } else {
        std::cerr << "Usage: " << argv[0] << " [-t|-m] <input_file>"
                  << std::endl;
        return -1;
    }
    /*** END PARAMETER PARSING ***/

    int nTh = -1;
    int nPh = -1;

    LegacyMHDConversion::getGridSize(inputFile, &nTh, &nPh, comm);
    if (nTh <= 0 || nPh <= 0) {
        throw std::runtime_error(
            "Error reading source term data, ending execution");
    }

    // p1 p2 p3
    // map
    // [1,2,3,4,5]
    auto map = rcp(new Map(nTh * nPh, 0, comm));
    auto sourceTerm = rcp(new Vector(map));
    RCP<Coordinates> coords;
    RCP<SolarModel> solarModel =
        rcp<SolarModel>(new SolarModel(YEAR, DAY, MONTH, HOUR));
    RCP<DipoleModel> dipoleModel = rcp<DipoleModel>(new DipoleModel());

    /*** BEGIN MHD INTERPOLATION ***/
    LegacyMHDConversion::processLegacyOutput(coords, dipoleModel, solarModel,
                                             sourceTerm, map, comm, nTh, nPh,
                                             inputFile);
    /*** END MHD INTERPOLATION ***/

    /*** BEGIN CONDUCTANCE SOLVE ***/

    auto [auroralConductance, euvConductance, conductance] =
        EmpConductance(coords, map, SIG0).computeConductance(KP, F107);
    /*** END CONDUCTANCE SOLVE ***/

    /*** BEGIN POTENTIAL SOLVE ***/
    TlSolver solver(coords, conductance, sourceTerm, map);
    VectorRCP result = solver.calculatePotential();
    /*** END POTENTIAL SOLVE ***/

    /*** BEGIN PLOTTING ***/
    if (useTecplot) {
        exportToTecplot("data/solvedPotential.dat", coords, result, sourceTerm,
                        auroralConductance, euvConductance, conductance, comm,
                        coords->nTh, coords->nPh);
    } else {
        exportToMatplotlib("data/solvedPotential.txt", coords->multiVector(),
                           result, comm, coords->nTh, coords->nPh);
    }
    return 0;
    /*** END PLOTTING ***/

    // TODO: Proper post-processing, including J and others
}
