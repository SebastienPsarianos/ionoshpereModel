#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/postProcessing/EField.hpp"
#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/tecplot.hpp"
#include "ionosphere/utils/Coordinates.hpp"
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
    double dTh = -1;
    double dPh = -1;

    LegacyMHDConversion::getGridSize(inputFile, &nTh, &nPh, comm);
    if (nTh <= 0 || nPh <= 0) {
        throw std::runtime_error(
            "Error reading source term data, ending execution");
    }

    auto map = rcp(new Map(nTh * nPh, 0, comm));
    auto coordinates = rcp(new MultiVector(map, 2));
    auto conductance = rcp(new MultiVector(map, 3));
    auto sourceTerm = rcp(new Vector(map));

    /*** BEGIN MHD INTERPOLATION ***/
    LegacyMHDConversion::processLegacyOutput(
        inputFile, dTh, dPh, coordinates, sourceTerm, comm, nTh, nPh, THETA0);
    /*** END MHD INTERPOLATION ***/

    /*** BEGIN CONDUCTANCE SOLVE ***/
    EmpConductance(conductance, coordinates, map, SIG0, nTh, nPh, YEAR, DAY,
                   MONTH, HOUR)
        .computeConductance(KP, F107);
    /*** END CONDUCTANCE SOLVE ***/

    /*** BEGIN POTENTIAL SOLVE ***/
    TlSolver test = TlSolver(nTh, nPh, dTh, dPh, map);
    VectorRCP result =
        test.calculatePotential(conductance, coordinates, sourceTerm);
    /*** END POTENTIAL SOLVE ***/

    /*** BEGIN PLOTTING ***/
    if (useTecplot) {
        exportToTecplot("../data/test3.dat", coordinates, result, sourceTerm,
                        comm, nTh, nPh);
    } else {
        exportToMatplotlib("../data/solvedPotential.txt", coordinates, result,
                           sourceTerm, comm, nTh, nPh);
    }
    return 0;
    /*** END PLOTTING ***/

    // TODO: Proper post-processing, including J and others
}
