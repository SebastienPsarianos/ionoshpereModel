#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/tecplot.hpp"
#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <Tpetra_Core.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Import.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_TpetraCrsGraphAdapter.hpp>
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

    /*** BEGIN ZOLTAN2 REPARTITIONING ***/
    // Build matrix on the naive map so Zoltan2 can analyze the graph
    TlSolver naiveSolver(coords, conductance, sourceTerm, map);
    auto naiveMatrix = naiveSolver.buildMatrix();

    // Use Zoltan2 hypergraph partitioning on the matrix sparsity pattern
    using GraphAdapter = Zoltan2::TpetraCrsGraphAdapter<
        Tpetra::CrsGraph<LocalOrd, GlobalOrd>>;

    auto graph = naiveMatrix->getCrsGraph();
    GraphAdapter adapter(rcp(new Tpetra::CrsGraph<LocalOrd, GlobalOrd>(*graph)));

    Teuchos::ParameterList zoltanParams;
    zoltanParams.set("algorithm", "phg");
    zoltanParams.set("num_global_parts", static_cast<int>(comm->getSize()));

    Zoltan2::PartitioningProblem<GraphAdapter> problem(&adapter, &zoltanParams);
    problem.solve();

    // Build new map from the Zoltan2 solution using Distributor for all-to-all
    auto& solution = problem.getSolution();
    auto partList = solution.getPartListView();
    auto myGids = map->getMyGlobalIndices();
    auto numLocal = map->getLocalNumElements();

    // Each rank sends its GIDs to the rank Zoltan2 assigned them to
    Teuchos::Array<int> exportPIDs(numLocal);
    Teuchos::Array<GlobalOrd> exportGIDs(numLocal);
    for (size_t i = 0; i < numLocal; i++) {
        exportPIDs[i] = static_cast<int>(partList[i]);
        exportGIDs[i] = myGids[i];
    }

    Tpetra::Distributor distributor(comm);
    size_t numImports = distributor.createFromSends(exportPIDs());

    Teuchos::Array<GlobalOrd> importGIDs(numImports);
    distributor.doPostsAndWaits(exportGIDs().getConst(), 1, importGIDs());

    auto invalidGST = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    auto newMap = rcp(new Map(invalidGST, importGIDs(), 0, comm));

    // Redistribute all data to the optimal map
    using Import = Tpetra::Import<LocalOrd, GlobalOrd>;
    Import importer(map, newMap);

    auto newCoordsMv = rcp(new MultiVector(newMap, 2));
    newCoordsMv->doImport(*coords->multiVector(), importer, Tpetra::INSERT);

    auto newConductance = rcp(new MultiVector(newMap, conductance->getNumVectors()));
    newConductance->doImport(*conductance, importer, Tpetra::INSERT);

    auto newSourceTerm = rcp(new Vector(newMap));
    newSourceTerm->doImport(*sourceTerm, importer, Tpetra::INSERT);

    auto newCoords = rcp(new Coordinates(newCoordsMv, dipoleModel, solarModel,
                                         nTh, nPh, coords->dTh, coords->dPh));
    /*** END ZOLTAN2 REPARTITIONING ***/

    /*** BEGIN POTENTIAL SOLVE ***/
    TlSolver solver(newCoords, newConductance, newSourceTerm, newMap);
    VectorRCP result = solver.calculatePotential();

    // Import result back to original map for output
    Import exportImporter(newMap, map);
    auto resultOnOrigMap = rcp(new Vector(map));
    resultOnOrigMap->doImport(*result, exportImporter, Tpetra::INSERT);
    result = resultOnOrigMap;
    /*** END POTENTIAL SOLVE ***/

    /*** BEGIN PLOTTING ***/
    if (useTecplot) {
        exportToTecplot("data/solvedPotential.dat", coords->multiVector(),
                        result, sourceTerm, auroralConductance, euvConductance,
                        conductance, comm, coords->nTh, coords->nPh);
    } else {
        exportToMatplotlib("data/solvedPotential.txt", coords->multiVector(),
                           result, comm, coords->nTh, coords->nPh);
    }
    return 0;
    /*** END PLOTTING ***/

    // TODO: Proper post-processing, including J and others
}
