#pragma once

#include <cmath>

#include "Teuchos_ParameterList.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_Comm.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

// TODO: move this to a cpp file
namespace Ionosphere {
inline void exportToTecplot(const Teuchos::ParameterList& ioParams,
                            Teuchos::RCP<Coordinates> coordWrapper,
                            Ionosphere::VectorRCP potential,
                            Ionosphere::VectorRCP sourceTerm,
                            Ionosphere::MultiVectorRCP auroralCondctance,
                            Ionosphere::MultiVectorRCP euvConductance,
                            Ionosphere::MultiVectorRCP conductance,
                            Teuchos::RCP<const Teuchos::Comm<int>> comm,
                            int nTh, int nPh) {

    std::string filename = ioParams.get<std::string>("output");

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
        // slower (phi). We write nPh+1 J-rows: nPh solver rows plus a closing
        // row that repeats ph=0, so Tecplot sees a closed surface with no gap
        // at the phi=0/2pi seam. The solver grid is unaffected.
        outFile << "ZONE I=" << nTh << ", J=" << nPh + 1
                << ", DATAPACKING=POINT\n";

        // --- Data Writing ---
        outFile << std::scientific << std::setprecision(6);

        auto writeRow = [&](int ph_src) {
            for (int th = 0; th < nTh; th++) {
                long long id = ph_src * nTh + th;
                double th_ = thVals[id], ph_ = phVals[id];
                double x = std::sin(th_) * std::cos(ph_);
                double y = std::sin(th_) * std::sin(ph_);
                double z = std::cos(th_);
                outFile << x << " " << y << " " << z << " " << th_ << " "
                        << ph_ << " " << potVals[id] << " " << jrVals[id]
                        << " " << pedersonVals[id] << " " << hallVals[id]
                        << " " << parallelVals[id] << " "
                        << auroralpedersonVals[id] << " "
                        << auroralhallVals[id] << " "
                        << auroralparallelVals[id] << " "
                        << euvpedersonVals[id] << " " << euvhallVals[id]
                        << " " << euvparallelVals[id] << " " << mltVals[id]
                        << "\n";
            }
        };

        for (int ph = 0; ph < nPh; ph++)
            writeRow(ph);

        // Closing row: repeat ph=0 to connect the phi=2pi seam in Tecplot
        writeRow(0);

        outFile.close();
        std::cout << "Successfully exported Tecplot data to: " << filename
                  << std::endl;
    }
}
} // namespace Ionosphere
