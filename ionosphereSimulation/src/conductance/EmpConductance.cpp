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

        // TODO: Testing
        double mlatDeg = mlat * 180.0 / M_PI;
        if (std::abs(mlatDeg) < 50.0) {
            pedersonConductances[i] = 0.0;
            hallConductances[i] = 0.0;
            continue;
        }

        double pedH0 = fourierSeries(pedersonData["max_latitude"], mlt);
        double hallH0 = fourierSeries(hallData["max_latitude"], mlt);

        double pedVal = epsteinFunction(
            std::abs(mlatDeg), fourierSeries(pedersonData["max_value"], mlt),
            pedH0, fourierSeries(pedersonData["up_slope"], mlt),
            fourierSeries(pedersonData["down_slope"], mlt));

        double hallVal = epsteinFunction(
            std::abs(mlatDeg), fourierSeries(hallData["max_value"], mlt),
            hallH0, fourierSeries(hallData["up_slope"], mlt),
            fourierSeries(hallData["down_slope"], mlt));

        // TODO: Testing
        if (std::abs(mlatDeg) < pedH0) {
            pedersonConductances[i] = std::max(0.0, pedVal);
        } else {
            pedersonConductances[i] = std::max(0.55, pedVal);
        }

        if (std::abs(mlatDeg) < hallH0) {
            hallConductances[i] = std::max(0.0, hallVal);
        } else {
            hallConductances[i] = std::max(0.55, hallVal);
        }
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

        // TODO: Testing
        hppPedersonConductance[i] = std::max(hppPedersonConductance[i], 0.25);
        hppHallConductance[i] = std::max(hppHallConductance[i], 0.25);

        hppParallelConductance[i] = _sig0;
    }
}
