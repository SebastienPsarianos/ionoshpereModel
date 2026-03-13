#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_Comm.hpp>
#include <Tpetra_Core.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace Ionosphere;
using Teuchos::rcp;
using Teuchos::RCP;

void LegacyMHDConversion::processLegacyOutput(RCP<Coordinates>& coordinates,
                                              RCP<DipoleModel> dipoleModel,
                                              RCP<SolarModel> solarModel,
                                              VectorRCP& sourceTerm, MapRCP map,
                                              CommRCP comm, int nTh, int nPh,
                                              const std::string& filename) {
    // Create root map so we only hit the file once
    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        nTh * nPh, (comm->getRank() == 0 ? nTh * nPh : 0), 0, comm));

    // Build the vectors for the root ran
    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootSourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));

    auto thVals = rootCoords->getDataNonConst(0);
    auto phVals = rootCoords->getDataNonConst(1);
    auto sourceVals = rootSourceTerm->getDataNonConst();

    double dTh = 0.0;
    double dPh = 0.0;

    if (comm->getRank() == 0) {
        std::fstream jrData = std::fstream(filename);
        if (!jrData.is_open()) {
            throw std::runtime_error("Failed opening radial current data");
        }

        std::string line;
        std::getline(jrData, line); // Ignore the line with the grid sizes

        for (int th = 0; th < nTh; th++) {
            for (int ph = 0; ph < nPh; ph++) {
                long long globalId = ph * nTh + th;
                jrData >> thVals[globalId] >> phVals[globalId] >>
                    sourceVals[globalId];

                if (!jrData) {
                    throw std::runtime_error(
                        "File structure incorrect for JR input data");
                }
            }
        }
        dTh = thVals[1] - thVals[0];
        dPh = phVals[nTh] - phVals[0];
    }

    auto coordVector = rcp(new Ionosphere::MultiVector(map, 2));
    sourceTerm = rcp(new Ionosphere::Vector(map));
    Tpetra::Export<int, long long> exporter(rootMap, map);

    coordVector->doExport(*rootCoords, exporter, Tpetra::INSERT);
    sourceTerm->doExport(*rootSourceTerm, exporter, Tpetra::INSERT);

    comm->broadcast(0, sizeof(double), (char*)&dTh);
    comm->broadcast(0, sizeof(double), (char*)&dPh);

    coordinates = rcp(new Coordinates(coordVector, dipoleModel, solarModel, nTh,
                                      nPh, dTh, dPh));
}
void LegacyMHDConversion::getGridSize(
    std::string filename, int* nTh, int* nPh,
    Teuchos::RCP<const Teuchos::Comm<int>> comm) {
    if (comm->getRank() == 0) {
        std::fstream dataFile = std::fstream(filename);
        if (!dataFile.is_open()) {
            throw std::runtime_error("Failed opening radial current data");
        }

        std::string line;
        if (std::getline(dataFile, line)) {
            std::sscanf(line.c_str(), "nTh: %d, nPh: %d", nTh, nPh);
        }

        std::cout << *nTh << " " << *nPh << std::endl;
    }

    comm->broadcast(0, sizeof(int), (char*)nPh);
    comm->broadcast(0, sizeof(int), (char*)nTh);
}
