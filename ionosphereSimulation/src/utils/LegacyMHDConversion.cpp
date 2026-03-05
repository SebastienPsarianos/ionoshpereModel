#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include "Kokkos_TeuchosCommAdapters.hpp"
#include "Teuchos_Comm.hpp"
#include <Tpetra_Core.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

void LegacyMHDConversion::processLegacyOutput(
    std::string filename, double& dTh, double& dPh,
    Teuchos::RCP<Tpetra::MultiVector<double, int, long long>> coords,
    Teuchos::RCP<Tpetra::Vector<double, int, long long>> sourceTerm,
    Teuchos::RCP<const Teuchos::Comm<int>> comm, int nTh, int nPh,
    double THETA0) {

    // Create root map so we only hit the file once
    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        nTh * nPh, (comm->getRank() == 0 ? nTh * nPh : 0), 0, comm));
    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootSourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));

    auto thVals = rootCoords->getDataNonConst(0);
    auto phVals = rootCoords->getDataNonConst(1);
    auto sourceVals = rootSourceTerm->getDataNonConst();

    if (comm->getRank() == 0) {
        std::fstream jrData = std::fstream(filename);
        if (!jrData.is_open()) {
            throw std::runtime_error("Failed opening radial current data");
        }

        std::string line;
        std::getline(jrData, line); // Ignore the line with the grid sizes

        for (int th = 0; th < nTh; th++) {
            for (int ph = 0; ph < nPh; ph++) {
                // if (ph == nPh) {
                //     double discard;
                //     jrData >> discard >> discard >> discard;
                //     // We don't store the last redundant ph value
                //     continue;
                // }

                long long globalId = ph * nTh + th;
                jrData >> thVals[globalId] >> phVals[globalId] >>
                    sourceVals[globalId];

                if (!jrData) {
                    throw std::runtime_error(
                        "File structure corrupted for JR input data");
                }

                if (th == 0) {
                    thVals[globalId] = thVals[globalId] + THETA0;
                } else if (th == nTh - 1) {
                    thVals[globalId] = thVals[globalId] - THETA0;
                }
            }
        }
        dTh = thVals[1] - (thVals[0] - THETA0);
        dPh = phVals[nTh] - phVals[0];
    }

    auto distMap = coords->getMap();
    Tpetra::Export<int, long long> exporter(rootMap, distMap);

    // Export the new Vectors
    coords->doExport(*rootCoords, exporter, Tpetra::INSERT);
    sourceTerm->doExport(*rootSourceTerm, exporter, Tpetra::INSERT);

    // Broadcast dTh, dPh values
    comm->broadcast(0, sizeof(double), (char*)&dPh);
    comm->broadcast(0, sizeof(double), (char*)&dTh);
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

        std::cout << nTh << " " << nPh << std::endl;
    }

    comm->broadcast(0, sizeof(int), (char*)nPh);
    comm->broadcast(0, sizeof(int), (char*)nTh);
}
