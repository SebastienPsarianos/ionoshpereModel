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
double HOUR = 5;

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
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return -1;
    }

    std::string inputFile = argv[2];
    /*** END PARAMETER PARSING ***/

    int nTh = -1;
    int nPh = -1;

    LegacyMHDConversion::getGridSize(inputFile, &nTh, &nPh, comm);
    if (nTh <= 0 || nPh <= 0) {
        throw std::runtime_error(
            "Error reading source term data, ending execution");
    }

    auto map = rcp(new Map(nTh * nPh, 0, comm));

    /*** BEGIN MHD INTERPOLATION ***/
    auto sourceTerm = rcp(new Vector(map));
    RCP<Coordinates> coords;
    RCP<SolarModel> solarModel =
        rcp<SolarModel>(new SolarModel(YEAR, DAY, MONTH, HOUR));
    RCP<DipoleModel> dipoleModel =
        rcp<DipoleModel>(new DipoleModel(comm, YEAR, DAY, MONTH, HOUR));

    LegacyMHDConversion::processLegacyOutput(coords, dipoleModel, solarModel,
                                             sourceTerm, map, comm, nTh, nPh,
                                             inputFile);
    /*** END MHD INTERPOLATION ***/

    /*** BEGIN CONDUCTANCE CALC***/
    auto [auroralConductance, euvConductance, conductance] =
        EmpConductance(coords, map, SIG0).computeConductance(KP, F107);
    /*** END CONDUCTANCE CALC ***/

    /*** BEGIN POTENTIAL SOLVE ***/
    TlSolver solver(coords, conductance, sourceTerm, map);
    VectorRCP result = solver.calculatePotential();
    /*** END POTENTIAL SOLVE ***/

    /*** BEGIN PLOTTING ***/
    exportToTecplot("data/solvedPotential.dat", coords, result, sourceTerm,
                    auroralConductance, euvConductance, conductance, comm,
                    coords->nTh, coords->nPh);

    return 0;
    /*** END PLOTTING ***/

    // TODO: Proper post-processing, including J and others
}
