#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <Tpetra_Core.hpp>
#include <fstream>
#include <stdexcept>

// TODO:
// - Move the grid by 0.5dTh, deleting the bottom point
// - Calculate dTh dPh
// - delete the last point in phi
void LegacyMHDConversion::processLegacyOutput(
    std::string filename, double& dTh, double& dPh,
    Teuchos::RCP<Tpetra::MultiVector<double, int, long long>> coords,
    Teuchos::RCP<Tpetra::Vector<double, int, long long>> sourceTerm,
    Teuchos::RCP<const Teuchos::Comm<int>> comm, int nTh, int nPh) {
    std::fstream jrData = std::fstream(filename);

    // Create root map so we only have one core reading the file
    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        nTh * nPh, (comm->getRank() == 0 ? nTh * nPh : 0), 0, comm));
    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootSourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));

    if (!jrData.is_open()) {
        throw std::runtime_error("Failed opening radial current data");
    }

    if (comm->getRank() == 0) {
        std::string line;
        std::getline(jrData, line); // Ignore the line with the grid sizes
        for (int th = 0; th < nTh; th++) {
            for (int ph = 0; ph < nPh; ph++) {
                if (std::getline(jrData, line)) {
                    double thetaVal, phiVal, radVal;

                    if (std::sscanf(line.c_str(), "%lf %lf %lf", &thetaVal,
                                    &phiVal, &radVal) != 3) {
                        throw std::runtime_error("Error parsing line" + line);
                    }

                    long long globalId = ph * nTh + th;

                    rootCoords->replaceGlobalValue(globalId, 0, thetaVal);
                    rootCoords->replaceGlobalValue(globalId, 1, phiVal);
                    rootSourceTerm->replaceGlobalValue(globalId, radVal);

                } else {
                    throw std::length_error(
                        "File length doesn't match provided coordinates");
                }
            }
        }
        dTh = rootCoords->getData(0)[1] - rootCoords->getData(0)[0];
        dPh = rootCoords->getData(1)[nPh] - rootCoords->getData(1)[0];
    }

    comm->broadcast(0, sizeof(int), (char*)&dPh);
    comm->broadcast(0, sizeof(int), (char*)&dTh);

    auto distMap = coords->getMap();
    Tpetra::Export<int, long long> exporter(rootMap, distMap);

    coords->doExport(*rootCoords, exporter, Tpetra::INSERT);
    sourceTerm->doExport(*rootSourceTerm, exporter, Tpetra::INSERT);
}

// TODO: Maybe move the broadcast into here
void LegacyMHDConversion::getGridSize(std::string filename, int* nTh,
                                      int* nPh) {
    std::fstream dataFile = std::fstream(filename);
    if (!dataFile.is_open()) {
        throw std::runtime_error("Failed opening radial current data");
    }

    std::string line;
    if (std::getline(dataFile, line)) {
        std::sscanf(line.c_str(), "nTh: %d, nPh: %d", nTh, nPh);
    }
}
