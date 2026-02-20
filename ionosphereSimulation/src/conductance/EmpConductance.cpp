#include "ionosphere/conductance/EmpConductance.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "ionosphere/conductance/utils.hpp"
#include "ionosphere/utils/Grid.hpp"

#include <cmath>
#include <fstream>
#include <nlohmann/json>
#include <string>

EmpConductance::EmpConductance(MultiVectorRcp& sigma, MultiVectorRcp& coords,
                               size_t nTh, size_t nPh, double sig0, int year,
                               int day, int month, int hour, MapRcp map)
    : _map(map), _sigma(sigma), _coords(coords),
      _euvConductance(new Tpetra::MultiVector<double, int, long long>(map, 3)),
      _auroralConductance(
          new Tpetra::MultiVector<double, int, long long>(map, 3)),
      _nTh(nTh), _nPh(nPh), _sig0(sig0) {

    _readAndSyncJson(map->getComm());
    computeDipoleRotationMatrix(_rotationMatrix);
    computeGrenaTimescales(_utTime, _ttTime, year, month, day, hour);
}

void EmpConductance::computeConductance(int kp, double f107) {

    // TODO: Set kp and f107 value dynamically
    _computeAuroralConductance(kp);
    _computeEuvConductance(f107);
    _computeHppConductance();
}

void EmpConductance::_readAndSyncJson(
    Teuchos::RCP<const Teuchos::Comm<int>> comm) {
    using nlohmann::json;

    const int myRank = comm->getRank();
    const int rootRank = 0;
    std::string serializedData;
    int dataLength = 0;

    if (myRank == rootRank) {
        try {
            std::ifstream pedFile(
                "../src/conductance/data/hardyPedersonCoeff.json");
            std::ifstream hallFile(
                "../src/conductance/data/hardyHallCoeff.json");

            if (!pedFile.is_open() || !hallFile.is_open()) {
                throw std::runtime_error(
                    "Unable to open hardy coefficient files");
            }

            json pedersonJson, hallJson;
            hallFile >> hallJson;
            pedFile >> pedersonJson;

            _coefficientJson["pederson"] = pedersonJson;
            _coefficientJson["hall"] = hallJson;

            serializedData = _coefficientJson.dump();
            dataLength = static_cast<int>(serializedData.size());

        } catch (const std::exception& e) {
            std::cerr << "JSON parse error" << e.what() << std::endl;
            dataLength = -1;
        }
    }

    Teuchos::broadcast(*comm, rootRank, 1, &dataLength);

    if (dataLength < 0)
        return;

    if (myRank != rootRank) {
        serializedData.resize(dataLength);
    }

    Teuchos::broadcast(*comm, rootRank, dataLength, &serializedData[0]);
    try {
        _coefficientJson = json::parse(serializedData);
    } catch (const std::exception& e) {
        if (myRank == rootRank) {
            std::cerr << "JSON parse error: " << e.what() << std::endl;
        }
    }
}

void EmpConductance::_computeEuvConductance(double f107) {
    using std::cos, std::sin, std::sqrt;
    auto pedersonConductances = _euvConductance->getDataNonConst(0);
    auto hallConductances = _euvConductance->getDataNonConst(1);
    auto thVals = _coords->getDataNonConst(0);
    auto phVals = _coords->getDataNonConst(1);

    for (int i = 0; i < thVals.size(); i++) {
        GeoSph currCoord = {.theta = thVals[i], .phi = phVals[i]};
        double sza = computeSolarZenith(_utTime, _ttTime, convert(currCoord));

        if (sza >= M_PI / 2) {
            pedersonConductances[i] = 0;
            hallConductances[i] = 0;
        } else {
            pedersonConductances[i] =
                pow(f107, 0.49) * (0.34 * cos(sza) + 0.93 * sqrt(cos(sza)));
            hallConductances[i] =
                pow(f107, 0.53) * (0.81 * cos(sza) + 0.54 * sqrt(cos(sza)));
        }
    }
}

// TODO: Figure out how to input the values equatorward of the max etc...
void EmpConductance::_computeAuroralConductance(int kp) {
    using nlohmann::json;
    if (kp > 6 || kp < 0)
        throw std::runtime_error("Invalid kp value supplied for auroral "
                                 "conductance, should be [0,6]");

    json pedersonData = _coefficientJson["pederson"]["K" + std::to_string(kp)];
    json hallData = _coefficientJson["hall"]["K" + std::to_string(kp)];

    auto pedersonConductances = _auroralConductance->getDataNonConst(0);
    auto hallConductances = _auroralConductance->getDataNonConst(1);
    auto thVals = _coords->getDataNonConst(0);
    auto phVals = _coords->getDataNonConst(1);

    for (int i = 0; i < thVals.size(); i++) {
        GeoSph currCoord = {.theta = thVals[i], .phi = phVals[i]};

        MagSph observerPositionMag;
        geoCentricToDipole(observerPositionMag, currCoord, _rotationMatrix);

        GeoGeo subsolarPosition;
        computeSubSolar(subsolarPosition, _utTime, _ttTime);
        MagSph subsolarPositionMag;
        geoCentricToDipole(subsolarPositionMag, convert(subsolarPosition),
                           _rotationMatrix);

        double mlt = computeMagneticLocalTime(convert(subsolarPositionMag),
                                              convert(observerPositionMag));

        double mlat = convert(observerPositionMag).latitude;

        pedersonConductances[i] =
            epsteinFunction(mlat, fourierSeries(pedersonData["max_value"], mlt),
                            fourierSeries(pedersonData["max_latitude"], mlt),
                            fourierSeries(pedersonData["up_slope"], mlt),
                            fourierSeries(pedersonData["down_slope"], mlt));
        hallConductances[i] =
            epsteinFunction(convert(observerPositionMag).latitude,
                            fourierSeries(hallData["max_value"], mlt),
                            fourierSeries(hallData["max_latitude"], mlt),
                            fourierSeries(hallData["up_slope"], mlt),
                            fourierSeries(hallData["down_slope"], mlt));
    }
}

void EmpConductance::_computeHppConductance() {
    auto pedersonConductancesAur = _auroralConductance->getDataNonConst(0);
    auto hallConductancesAur = _auroralConductance->getDataNonConst(1);
    auto pedersonConductancesEuv = _euvConductance->getDataNonConst(0);
    auto hallConductancesEuv = _euvConductance->getDataNonConst(1);
    auto thVals = _coords->getDataNonConst(0);
    auto phVals = _coords->getDataNonConst(1);

    auto hppPedersonConductance = _sigma->getDataNonConst(0);
    auto hppHallConductance = _sigma->getDataNonConst(1);
    auto hppParallelConductance = _sigma->getDataNonConst(2);

    for (int i = 0; i < thVals.size(); i++) {
        hppHallConductance[i] =
            sqrt(hallConductancesAur[i] * hallConductancesAur[i] +
                 hallConductancesEuv[i] * hallConductancesEuv[i]);
        hppPedersonConductance[i] =
            sqrt(pedersonConductancesAur[i] * pedersonConductancesAur[i] +
                 pedersonConductancesEuv[i] * pedersonConductancesEuv[i]);

        hppParallelConductance[i] = _sig0;
    }
}
void EmpConductance::_calcSigmaDer(Grid<DSigma>& dConductance,
                                   Grid<Sigma>& sigma, Grid<GeoSph>& coords) {

    const double dTh = coords(1, 0).theta - coords(0, 0).theta;
    const double dPh = coords(0, 1).phi - coords(0, 0).phi;

    // TODO: Figure out pole behaviour for th
    for (size_t th = 1; th < _nTh - 1; th++) {
        // NOTE: Boundary points for phi should be continuous
        //      and wrap around to the start of the phi grid.
        //     Boundary points (phi = 0 and phi = _nPh - 1) are
        //      calculated as derivatives between phi = 1 and nPh-2
        dConductance(th, 0).dthph_ph =
            (sigma(th, 1).thph - sigma(th, _nPh - 2).thph) / (2 * dPh);
        dConductance(th, _nPh - 1).dthph_ph = dConductance(th, 0).dthph_ph;

        dConductance(th, 0).dphph_ph =
            (sigma(th, 1).phph - sigma(th, _nPh - 2).phph) / (2 * dPh);
        dConductance(th, _nPh - 1).dphph_ph = dConductance(th, 0).dphph_ph;

        // We also set the th derivatives to be equal accross the circular
        // boundary
        dConductance(th, 0).dthth_th =
            (sigma(th + 1, 0).thth - sigma(th - 1, 0).thth) / (2 * dTh);
        dConductance(th, _nPh - 1).dthth_th = dConductance(th, 0).dthth_th;

        dConductance(th, 0).dthph_th =
            (sigma(th + 1, 0).thph - sigma(th - 1, 0).thph) / (2 * dTh);
        dConductance(th, _nPh - 1).dthph_th = dConductance(th, 0).dthph_th;

        for (size_t ph = 1; ph < _nPh - 1; ph++) {
            dConductance(th, ph).dthth_th =
                (sigma(th + 1, ph).thth - sigma(th - 1, ph).thth) / (2 * dTh);

            dConductance(th, ph).dthph_ph =
                (sigma(th, ph + 1).thph - sigma(th, ph - 1).thph) / (2 * dPh);

            dConductance(th, ph).dthph_th =
                (sigma(th + 1, ph).thph - sigma(th - 1, ph).thph) / (2 * dTh);

            dConductance(th, ph).dphph_ph =
                (sigma(th, ph + 1).phph - sigma(th, ph - 1).phph) / (2 * dPh);
        }
    }
}
