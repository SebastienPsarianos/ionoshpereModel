#include <cmath>

#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_Comm.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <fstream>
#include <iomanip>

inline void exportToTecplot(const std::string& filename,
                            Teuchos::RCP<Coordinates> coordWrapper,
                            Ionosphere::VectorRCP potential,
                            Ionosphere::VectorRCP sourceTerm,
                            Ionosphere::MultiVectorRCP auroralCondctance,
                            Ionosphere::MultiVectorRCP euvConductance,
                            Ionosphere::MultiVectorRCP conductance,
                            Teuchos::RCP<const Teuchos::Comm<int>> comm,
                            int nTh, int nPh) {

    // 1. Create a root map where Rank 0 owns all the elements
    long long globalElements = nTh * nPh;
    size_t localElements = (comm->getRank() == 0) ? globalElements : 0;
    auto coords = coordWrapper->multiVector();
    auto distMap = coords->getMap();

    // Compute MLT on the distributed map
    auto mltVec =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(distMap));
    {
        auto mltData = mltVec->getDataNonConst(0);
        for (size_t i = 0; i < mltVec->getLocalLength(); i++) {
            mltData[i] = coordWrapper->localIdx2Mlt(i);
        }
    }

    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        globalElements, localElements, 0, comm));

    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootPotential =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));
    auto rootSourceTerm =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));
    auto rootConductance = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 3));
    auto rootAuroralConductance = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 3));
    auto rootEuvConductance = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 3));
    auto rootMlt =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));
    Tpetra::Export<int, long long> exporter(distMap, rootMap);

    rootCoords->doExport(*coords, exporter, Tpetra::INSERT);
    rootPotential->doExport(*potential, exporter, Tpetra::INSERT);
    rootSourceTerm->doExport(*sourceTerm, exporter, Tpetra::INSERT);
    rootConductance->doExport(*conductance, exporter, Tpetra::INSERT);
    rootAuroralConductance->doExport(*auroralCondctance, exporter,
                                     Tpetra::INSERT);
    rootEuvConductance->doExport(*euvConductance, exporter, Tpetra::INSERT);
    rootMlt->doExport(*mltVec, exporter, Tpetra::INSERT);

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
        auto pedersonVals = rootConductance->getData(0);
        auto hallVals = rootConductance->getData(1);
        auto parallelVals = rootConductance->getData(2);

        auto auroralpedersonVals = rootAuroralConductance->getData(0);
        auto auroralhallVals = rootAuroralConductance->getData(1);
        auto auroralparallelVals = rootAuroralConductance->getData(2);

        auto euvpedersonVals = rootEuvConductance->getData(0);
        auto euvhallVals = rootEuvConductance->getData(1);
        auto euvparallelVals = rootEuvConductance->getData(2);

        auto mltVals = rootMlt->getData(0);

        // --- Tecplot Header ---
        outFile << "TITLE = \"Ionosphere Potential Solution\"\n";
        outFile << "VARIABLES = \"X\", \"Y\", \"Z\", \"Theta\", \"Phi\", "
                   "\"Potential\", \"J_R\", \"Sigma_P\", \"Sigma_H\", "
                   "\"Sigma_0\", \"Sigma_P_AUR\", \"Sigma_H_AUR\", "
                   "\"Sigma_0_AUR\", "
                   "\"Sigma_P_EUV\", \"Sigma_H_EUV\", \"Sigma_0_EUV\", \"MLT\""
                   "\n";
        // Zone definition: I is the fastest varying index (theta), J is the
        // slower (phi)
        outFile << "ZONE I=" << nTh << ", J=" << nPh << ", DATAPACKING=POINT\n";

        // --- Data Writing ---
        outFile << std::scientific << std::setprecision(6);

        for (int ph = 0; ph < nPh; ph++) {
            for (int th = 0; th < nTh; th++) {
                long long id = ph * nTh + th;
                double th_ = thVals[id], ph_ = phVals[id];
                double x = std::sin(th_) * std::cos(ph_);
                double y = std::sin(th_) * std::sin(ph_);
                double z = std::cos(th_);
                outFile << x << " " << y << " " << z << " " << th_ << " " << ph_
                        << " " << potVals[id] << " " << jrVals[id] << " "
                        << pedersonVals[id] << " " << hallVals[id] << " "
                        << parallelVals[id] << " " << auroralpedersonVals[id]
                        << " " << auroralhallVals[id] << " "
                        << auroralparallelVals[id] << " " << euvpedersonVals[id]
                        << " " << euvhallVals[id] << " " << euvparallelVals[id]
                        << " " << mltVals[id] << "\n";
            }
        }

        outFile.close();
        std::cout << "Successfully exported Tecplot data to: " << filename
                  << std::endl;
    }
}

inline void exportToMatplotlib(
    const std::string& filename,
    Teuchos::RCP<const Tpetra::MultiVector<double, int, long long>> coords,
    Teuchos::RCP<const Tpetra::Vector<double, int, long long>> potential,
    Teuchos::RCP<const Teuchos::Comm<int>> comm, int nTh, int nPh) {

    long long globalElements = nTh * nPh;
    size_t localElements = (comm->getRank() == 0) ? globalElements : 0;

    auto rootMap = Teuchos::rcp(new Tpetra::Map<int, long long>(
        globalElements, localElements, 0, comm));

    auto rootCoords = Teuchos::rcp(
        new Tpetra::MultiVector<double, int, long long>(rootMap, 2));
    auto rootPotential =
        Teuchos::rcp(new Tpetra::Vector<double, int, long long>(rootMap));

    auto distMap = coords->getMap();
    Tpetra::Export<int, long long> exporter(distMap, rootMap);

    rootCoords->doExport(*coords, exporter, Tpetra::INSERT);
    rootPotential->doExport(*potential, exporter, Tpetra::INSERT);

    if (comm->getRank() == 0) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            throw std::runtime_error(
                "Failed to open file for matplotlib export!");
        }

        auto thVals = rootCoords->getData(0);
        auto phVals = rootCoords->getData(1);
        auto potVals = rootPotential->getData(0);

        outFile << std::scientific << std::setprecision(6);

        // Theta-outer, phi-inner order so plot2.py reshape((n_th, n_ph)) is
        // correct
        for (int th = 0; th < nTh; th++) {
            for (int ph = 0; ph < nPh; ph++) {
                long long id = ph * nTh + th;
                double th_ = thVals[id], ph_ = phVals[id];
                double x = std::sin(th_) * std::cos(ph_);
                double y = std::sin(th_) * std::sin(ph_);
                double z = std::cos(th_);
                outFile << x << " " << y << " " << z << " " << th_ << " " << ph_
                        << " " << potVals[id] << "\n";
            }
        }

        outFile.close();
        std::cout << "Successfully exported matplotlib data to: " << filename
                  << std::endl;
    }
}
