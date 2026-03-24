#include "ionosphere/solver/PoleDivThmBC.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ionosphere/TrilinosAliases.hpp"

using namespace Ionosphere;

PoleDivThmBC::PoleDivThmBC(const Teuchos::ParameterList& equationParams,
                           Teuchos::RCP<Coordinates> coords,
                           MultiVectorRCP conductance, MapRCP map)
    : _conductance(conductance), _map(map), _coords(coords),
      _radiusIonosphere2((equationParams.get<double>("ionosphere_altitude_m") +
                          equationParams.get<double>("earth_radius_m")) *
                         (equationParams.get<double>("ionosphere_altitude_m") +
                          equationParams.get<double>("earth_radius_m"))) {}

void PoleDivThmBC::apply(Ionosphere::MatrixRCP matrix,
                         Ionosphere::VectorRCP rhs) {
    _applyToMatrix(matrix);
    _applyToRHS(rhs);
}

// TODO: Lot of verfication to do here
void PoleDivThmBC::_applyToMatrix(Ionosphere::MatrixRCP matrix) {
    using std::sin;

    auto pedersonVals = _conductance->getDataNonConst(0);
    auto hallVals = _conductance->getDataNonConst(1);
    auto parallelVals = _conductance->getDataNonConst(2);

    auto coordsMv = _coords->multiVector();
    auto th = coordsMv->getDataNonConst(0);
    auto ph = coordsMv->getDataNonConst(1);

    auto myGridPoints = _map->getMyGlobalIndices();

    auto poleData = _gatherPoleData();
    auto sigThTh = poleData->getData(1);
    auto sigThPh = poleData->getData(2);

    // Set phi == 0 on the equator to be a gauge condition
    for (GlobalOrd i = 0; i < static_cast<GlobalOrd>(myGridPoints.size());
         i++) {
        GlobalOrd gridPoint = myGridPoints[i];
        auto [theta, phi] = _coords->globalIdx2ThetaPhi(gridPoint);

        Teuchos::Array<GlobalOrd> matrixIdcs;
        Teuchos::Array<Scalar> vals;

        // Pin the gauge point on the equator and then skip
        if (theta == _coords->nTh / 2 && phi == 0) {
            matrixIdcs.push_back(gridPoint);
            vals.push_back(1);
            matrix->insertGlobalValues(gridPoint, matrixIdcs, vals);
            continue;
        }

        GlobalOrd primaryPoleId = (theta == 0) ? 0 : (_coords->nTh - 1);

        // Applying the polar cap flux boundary condition (see notes)
        if (theta == 0 || theta == _coords->nTh - 1) {
            if (phi == 0) {
                Scalar theta0 = _coords->dTh / 2;

                // Grab the indices in the poleData multivector
                LocalOrd pole = theta == 0 ? 0 : (_coords->nPh + 1);
                LocalOrd ringStart = (theta == 0) ? 1 : (_coords->nPh + 2);

                // Multiply terms by negative 1 at the south pole to
                // account for  negative normal vector direction
                Scalar multiplier = theta == 0 ? 1 : -1;

                // Summing average conductance contributions across all phi
                // for pole cap (taking an average of the pole and the first
                // ring of theta values). Here we are using the 0 point for
                // consistency sake.
                Scalar poleCoefficient = 0;
                for (GlobalOrd j = 0; j < _coords->nPh; j++) {
                    poleCoefficient +=
                        (sigThTh[ringStart + j] + sigThTh[pole]) / 2;
                }
                poleCoefficient *= -sin(theta0) * _coords->dPh / _coords->dTh;

                matrixIdcs.push_back(gridPoint);
                vals.push_back(poleCoefficient);

                // Handle the contributions from the surrounding ring
                for (GlobalOrd j = 0; j < _coords->nPh; j++) {
                    Scalar sigThThAvg =
                        (sigThTh[ringStart + j] + sigThTh[pole]) / 2;

                    GlobalOrd prevIdx = (j - 1 + _coords->nPh) % _coords->nPh;
                    GlobalOrd nextIdx = (j + 1) % _coords->nPh;
                    Scalar sigThPhAvg = (sigThPh[ringStart + prevIdx] -
                                         sigThPh[ringStart + nextIdx]) /
                                        8;

                    Scalar ringPointCoeff =

                        (_coords->dPh * sin(theta0) / _coords->dTh) *
                            sigThThAvg +
                        multiplier * sigThPhAvg;

                    GlobalOrd ringTheta = (theta == 0) ? 1 : (_coords->nTh - 2);
                    GlobalOrd ringGlobIdx = j * _coords->nTh + ringTheta;

                    matrixIdcs.push_back(ringGlobIdx);
                    vals.push_back(ringPointCoeff);
                }

                matrix->insertGlobalValues(gridPoint, matrixIdcs, vals);

            } else {
                // Just constrain to phi == 0
                matrixIdcs.push_back(gridPoint);
                vals.push_back(1);
                matrixIdcs.push_back(primaryPoleId);
                vals.push_back(-1);
                matrix->insertGlobalValues(gridPoint, matrixIdcs, vals);
            }
        }
    }
}

// TODO: Verify
void PoleDivThmBC::_applyToRHS(VectorRCP rhs) {
    GlobalOrd pinPoint = _coords->nTh / 2;
    {
        auto rhsData = rhs->getDataNonConst();
        auto myGridPoints = _map->getMyGlobalIndices();
        for (LocalOrd i = 0; i < static_cast<LocalOrd>(myGridPoints.size());
             i++) {
            GlobalOrd currentGid = myGridPoints[i];
            GlobalOrd theta = currentGid % _coords->nTh;
            GlobalOrd phi = currentGid / _coords->nTh;

            if (currentGid == pinPoint) {
                rhsData[i] = 0.0;
            } else if (theta == 0 || theta == _coords->nTh - 1) {
                // We set the two primary pole points to the cap condition
                if (phi == 0) {
                    Scalar theta0 = _coords->dTh / 2.0;
                    Scalar jR_pole = rhsData[i] / _radiusIonosphere2;
                    rhsData[i] = jR_pole * _radiusIonosphere2 * 2.0 * M_PI *
                                 (1.0 - std::cos(theta0));
                } else {
                    // we constrain the rest of the pole points to theta
                    rhsData[i] = 0.0;
                }
            }
        }
    }
}

// TODO: Lot of verification to do here
MultiVectorRCP PoleDivThmBC::_gatherPoleData() {
    auto comm = _map->getComm();
    size_t localSize = 2 * (_coords->nPh + 1);

    // Create a local map that contains all pole values
    auto localMap = Teuchos::rcp(
        new Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                localSize, 0, comm));

    // 0 = theta, 1 = sigThTh, 2 = sigThPh (needed for the pole boundary
    // condition)
    auto poleData = Teuchos::rcp(new MultiVector(localMap, 3));
    poleData->putScalar(0.0);

    auto pdTh = poleData->getDataNonConst(0);
    auto sigThTh = poleData->getDataNonConst(1);
    auto sigThPh = poleData->getDataNonConst(2);

    auto coordsMv = _coords->multiVector();
    auto thVals = coordsMv->getDataNonConst(0);
    auto pedVals = _conductance->getDataNonConst(0);
    auto hallVals = _conductance->getDataNonConst(1);
    auto parVals = _conductance->getDataNonConst(2);

    // Construct lookup table of the entries we want to keep
    // 0 is the pole and then the rest are the ring nodes. pdIdx gives the index
    // in the localVector. TODO: This needs to be cleaner
    struct PoleEntry {
        GlobalOrd currentGid;
        GlobalOrd pdIdx;
    };

    std::vector<PoleEntry> entries;

    // Add the north pole
    entries.push_back({0, 0});
    // Add all of the north pole ring points
    for (GlobalOrd j = 0; j < _coords->nPh; j++)
        entries.push_back(
            {static_cast<GlobalOrd>(j * _coords->nTh + 1), j + 1});

    // Add the south pole
    entries.push_back({_coords->nTh - 1, _coords->nPh + 1});

    // Add the south pole ring points
    for (GlobalOrd j = 0; j < _coords->nPh; j++)
        entries.push_back(
            {j * _coords->nTh + _coords->nTh - 2, _coords->nPh + 2 + j});

    // Compute the required spherical conductance values for each pole point
    for (const auto& entry : entries) {
        if (!_map->isNodeGlobalElement(entry.currentGid))
            continue;

        auto li = _map->getLocalElement(entry.currentGid);
        Scalar th = thVals[li];
        Scalar cosE = (-2.0 * std::cos(th)) /
                      std::sqrt(1.0 + 3.0 * std::cos(th) * std::cos(th));
        Scalar sinE =
            std::sin(th) / std::sqrt(1.0 + 3.0 * std::cos(th) * std::cos(th));
        Scalar C = parVals[li] * cosE * cosE + pedVals[li] * sinE * sinE;

        pdTh[entry.pdIdx] = th;
        sigThTh[entry.pdIdx] = parVals[li] * pedVals[li] / C;
        sigThPh[entry.pdIdx] = parVals[li] * hallVals[li] * (-cosE) / C;
    }

    // Use a reduce to get the values available in all processors
    auto poleDataLocal = Teuchos::rcp(new MultiVector(*poleData));
    for (size_t col = 0; col < 3; col++) {
        auto src = poleDataLocal->getDataNonConst(col);
        auto dst = poleData->getDataNonConst(col);
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM,
                           static_cast<int>(localSize), src.getRawPtr(),
                           dst.getRawPtr());
    }
    return poleData;
}
