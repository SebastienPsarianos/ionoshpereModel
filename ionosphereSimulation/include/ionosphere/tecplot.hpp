#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_Vector_decl.hpp"
#include <Teuchos_Comm.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Export.hpp>

#include <fstream>

#include <iomanip>

void exportToTecplot(
    const std::string& filename,
    Teuchos::RCP<const Tpetra::MultiVector<double, int, long long>> coords,
    Teuchos::RCP<const Tpetra::Vector<double, int, long long>> potential,
    Teuchos::RCP<const Tpetra::Vector<double, int, long long>> sourceTerm,
    Teuchos::RCP<const Teuchos::Comm<int>> comm, int nTh, int nPh) {

    // 1. Create a root map where Rank 0 owns all the elements
    long long globalElements = nTh * nPh;
    size_t localElements = (comm->getRank() == 0) ? globalElements : 0;

    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        globalElements, localElements, 0, comm));

    // 2. Create the root vectors to hold the gathered data
    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootPotential =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));
    auto rootSourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));

    // 3. Export the distributed data to Rank 0
    auto distMap = coords->getMap();
    Tpetra::Export<int, long long> exporter(distMap, rootMap);

    rootCoords->doExport(*coords, exporter, Tpetra::INSERT);
    rootPotential->doExport(*potential, exporter, Tpetra::INSERT);
    rootSourceTerm->doExport(*sourceTerm, exporter, Tpetra::INSERT);

    // 4. Rank 0 writes the Tecplot file
    if (comm->getRank() == 0) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            throw std::runtime_error("Failed to open file for Tecplot export!");
        }

        // Extract the raw data arrays from the root vectors
        auto thVals = rootCoords->getData(0);
        auto phVals = rootCoords->getData(1);
        auto potVals = rootPotential->getData(0);
        auto jrVals = rootSourceTerm->getData(0);

        // --- Tecplot Header ---
        outFile << "TITLE = \"Ionosphere Potential Solution\"\n";
        outFile << "VARIABLES = \"Theta\", \"Phi\", \"Potential\", \"J_R\"\n";
        // Zone definition: I is the fastest varying index (theta), J is the
        // slower (phi)
        outFile << "ZONE I=" << nTh << ", J=" << nPh << ", DATAPACKING=POINT\n";

        // --- Data Writing ---
        outFile << std::scientific << std::setprecision(6);

        for (int ph = 0; ph < nPh; ph++) {
            for (int th = 0; th < nTh; th++) {
                long long id = ph * nTh + th;
                outFile << thVals[id] << " " << phVals[id] << " " << potVals[id]
                        << " " << jrVals[id] << "\n";
            }
        }

        outFile.close();
        std::cout << "Successfully exported Tecplot data to: " << filename
                  << std::endl;
    }
}

void exportToMatplotlib(
    const std::string& filename,
    Teuchos::RCP<const Tpetra::MultiVector<double, int, long long>> coords,
    Teuchos::RCP<const Tpetra::Vector<double, int, long long>> potential,
    Teuchos::RCP<const Tpetra::Vector<double, int, long long>> sourceTerm,
    Teuchos::RCP<const Teuchos::Comm<int>> comm, int nTh, int nPh) {

    long long globalElements = nTh * nPh;
    size_t localElements = (comm->getRank() == 0) ? globalElements : 0;

    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        globalElements, localElements, 0, comm));

    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootPotential =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));
    auto rootSourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));

    auto distMap = coords->getMap();
    Tpetra::Export<int, long long> exporter(distMap, rootMap);

    rootCoords->doExport(*coords, exporter, Tpetra::INSERT);
    rootPotential->doExport(*potential, exporter, Tpetra::INSERT);
    rootSourceTerm->doExport(*sourceTerm, exporter, Tpetra::INSERT);

    if (comm->getRank() == 0) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            throw std::runtime_error(
                "Failed to open file for matplotlib export!");
        }

        auto thVals = rootCoords->getData(0);
        auto phVals = rootCoords->getData(1);
        auto potVals = rootPotential->getData(0);
        auto jrVals = rootSourceTerm->getData(0);

        outFile << std::scientific << std::setprecision(6);

        // Theta-outer, phi-inner order so plot2.py reshape((n_th, n_ph)) is correct
        for (int th = 0; th < nTh; th++) {
            for (int ph = 0; ph < nPh; ph++) {
                long long id = ph * nTh + th;
                outFile << thVals[id] << " " << phVals[id] << " "
                        << potVals[id] << " " << jrVals[id] << "\n";
            }
        }

        outFile.close();
        std::cout << "Successfully exported matplotlib data to: " << filename
                  << std::endl;
    }
}
