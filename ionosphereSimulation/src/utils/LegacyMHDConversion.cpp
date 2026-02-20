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

    std::vector<double> thVals = std::vector<double>((nPh + 1) * (nTh + 1));
    std::vector<double> phVals = std::vector<double>((nPh + 1) * (nTh + 1));
    std::vector<double> sourceVals = std::vector<double>((nPh + 1) * (nTh + 1));

    if (comm->getRank() == 0) {
        std::string line;
        std::getline(jrData, line); // Ignore the line with the grid sizes
        for (int th = 0; th < nTh + 1; th++) {
            for (int ph = 0; ph < nPh + 1; ph++) {
                if (std::getline(jrData, line)) {

                    long long globalId = ph * nTh + th;
                    if (std::sscanf(line.c_str(), "%lf %lf %lf",
                                    &thVals[globalId], &phVals[globalId],
                                    &sourceVals[globalId]) != 3) {
                        throw std::runtime_error("Error parsing line" + line);
                    }
                } else {
                    throw std::length_error(
                        "File length doesn't match provided coordinates");
                }
            }
        }
        dTh = thVals[1] - thVals[0];
        dPh = phVals[nPh] - phVals[0];
    }

    for (int th = 0; th < (nTh + 1); th++) {
        for (int ph = 0; ph < nPh + 1; ph++) {
            long long globalId = ph * nTh + th;
        }
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
    int tempNTh = 0;
    int tempNPh = 0;

    std::string line;
    if (std::getline(dataFile, line)) {
        std::sscanf(line.c_str(), "nTh: %d, nPh: %d", &tempNTh, &tempNPh);
    }
    *nTh = tempNTh - 1;
    *nPh = tempNPh - 1;
}
