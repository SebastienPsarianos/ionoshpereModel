#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/conductance/utils.hpp"

#include <cmath>
#include <fstream>
#include <nlohmann/json>
#include <string>

using namespace Ionosphere;
using nlohmann::json;

EmpConductance::EmpConductance(Teuchos::RCP<Coordinates> coords, MapRCP map,
                               Scalar sig0)
    : _map(map), _euvConductance(new Ionosphere::MultiVector(map, 3)),
      _auroralConductance(new Ionosphere::MultiVector(map, 3)),
      _conductance(new Ionosphere::MultiVector(map, 3)), _sig0(sig0),
      _coords(coords) {

    _readAndSyncJson(map->getComm());
}

std::tuple<MultiVectorRCP, MultiVectorRCP, MultiVectorRCP>
EmpConductance::computeConductance(int kp, double f107) {
    _computeAuroralConductance(kp);
    _computeEuvConductance(f107);
    _computeHppConductance();

    return {_auroralConductance, _euvConductance, _conductance};
}

void EmpConductance::_readAndSyncJson(CommRCP comm) {

    const int myRank = comm->getRank();
    const int rootRank = 0;
    std::string serializedData;
    int dataLength = 0;

    if (myRank == rootRank) {
        try {
            std::ifstream pedFile(
                "src/conductance/data/hardyPedersonCoeff.json");
            std::ifstream hallFile("src/conductance/data/hardyHallCoeff.json");

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

    for (size_t i = 0; i < _coords->multiVector()->getLocalLength(); i++) {
        double sza = _coords->localIdx2Sza(i);

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

void EmpConductance::_computeAuroralConductance(int kp) {

    if (kp > 6 || kp < 0)
        throw std::runtime_error("Invalid kp value supplied for auroral "
                                 "conductance, should be [0,6]");

    json pedersonData = _coefficientJson["pederson"]["K" + std::to_string(kp)];
    json hallData = _coefficientJson["hall"]["K" + std::to_string(kp)];

    auto pedersonConductances = _auroralConductance->getDataNonConst(0);
    auto hallConductances = _auroralConductance->getDataNonConst(1);

    for (size_t i = 0; i < _coords->multiVector()->getLocalLength(); i++) {
        double mlt = _coords->localIdx2Mlt(i);
        // TODO: Make sure I should be using the abs here
        double mlat = std::abs(_coords->localIdx2MagGeo(i).latitude);

        // Clamping below 50 degrees latitude
        double mlatDeg = mlat * 180.0 / M_PI;
        if (mlatDeg < 50.0) {
            pedersonConductances[i] = 0.0;
            hallConductances[i] = 0.0;
            continue;
        }

        double pedH0 = fourierSeries(pedersonData["max_latitude"], mlt);
        double hallH0 = fourierSeries(hallData["max_latitude"], mlt);

        double pedVal = epsteinFunction(
            mlatDeg, fourierSeries(pedersonData["max_value"], mlt), pedH0,
            fourierSeries(pedersonData["up_slope"], mlt),
            fourierSeries(pedersonData["down_slope"], mlt));

        double hallVal =
            epsteinFunction(mlatDeg, fourierSeries(hallData["max_value"], mlt),
                            hallH0, fourierSeries(hallData["up_slope"], mlt),
                            fourierSeries(hallData["down_slope"], mlt));

        // The clamping described in the paper
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

    auto hppPedersonConductance = _conductance->getDataNonConst(0);
    auto hppHallConductance = _conductance->getDataNonConst(1);
    auto hppParallelConductance = _conductance->getDataNonConst(2);

    for (size_t i = 0; i < _coords->multiVector()->getLocalLength(); i++) {
        hppHallConductance[i] =
            sqrt(hallConductancesAur[i] * hallConductancesAur[i] +
                 hallConductancesEuv[i] * hallConductancesEuv[i]);
        hppPedersonConductance[i] =
            sqrt(pedersonConductancesAur[i] * pedersonConductancesAur[i] +
                 pedersonConductancesEuv[i] * pedersonConductancesEuv[i]);

        hppParallelConductance[i] = _sig0;

        // TODO: Testing minimum background conductance
        hppPedersonConductance[i] = std::max(hppPedersonConductance[i], 0.25);
        hppHallConductance[i] = std::max(hppHallConductance[i], 0.25);
        hppParallelConductance[i] = _sig0;
    }
}
