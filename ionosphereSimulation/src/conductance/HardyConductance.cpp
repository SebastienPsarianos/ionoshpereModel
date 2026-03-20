#include "ionosphere/conductance/HardyConductance.hpp"

#include <cmath>
#include <fstream>
#include <nlohmann/json>
#include <string>

using namespace Ionosphere;
using nlohmann::json;

HardyConductance::HardyConductance(Teuchos::RCP<Coordinates> coords, MapRCP map)
    : _map(map), _auroralConductance(new Ionosphere::MultiVector(map, 3)),
      _coords(coords) {

    _readAndSyncJson(map->getComm());
}

void HardyConductance::_readAndSyncJson(CommRCP comm) {
    const int myRank = comm->getRank();
    const int rootRank = 0;
    std::string serializedData;
    int dataLength = 0;

    if (myRank == rootRank) {
        try {
            std::ifstream pedFile("data/hardyPedersonCoeff.json");
            std::ifstream hallFile("data/hardyHallCoeff.json");

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

MultiVectorRCP HardyConductance::computeAuroralConductance(int kp) {

    if (kp > 6 || kp < 0)
        throw std::runtime_error("Invalid kp value supplied for auroral "
                                 "conductance, should be [0,6]");

    json pedersonData = _coefficientJson["pederson"]["K" + std::to_string(kp)];
    json hallData = _coefficientJson["hall"]["K" + std::to_string(kp)];

    auto pedersonConductances = _auroralConductance->getDataNonConst(0);
    auto hallConductances = _auroralConductance->getDataNonConst(1);

    for (size_t i = 0; i < _coords->multiVector()->getLocalLength(); i++) {
        double mlt = _coords->localIdx2Mlt(i);
        double mlat = std::abs(_coords->localIdx2MagGeo(i).latitude);

        // Clamping below 50 degrees latitude
        double mlatDeg = mlat * 180.0 / M_PI;
        if (mlatDeg < 50.0) {
            pedersonConductances[i] = 0.0;
            hallConductances[i] = 0.0;
            continue;
        }

        double pedH0 = _fourierSeries(pedersonData["max_latitude"], mlt);
        double hallH0 = _fourierSeries(hallData["max_latitude"], mlt);

        double pedVal = _epsteinFunction(
            mlatDeg, _fourierSeries(pedersonData["max_value"], mlt), pedH0,
            _fourierSeries(pedersonData["up_slope"], mlt),
            _fourierSeries(pedersonData["down_slope"], mlt));

        double hallVal = _epsteinFunction(
            mlatDeg, _fourierSeries(hallData["max_value"], mlt), hallH0,
            _fourierSeries(hallData["up_slope"], mlt),
            _fourierSeries(hallData["down_slope"], mlt));

        // The clamping described in the paper
        if (mlatDeg < pedH0) {
            pedersonConductances[i] = std::max(0.0, pedVal);
        } else {
            pedersonConductances[i] = std::max(0.55, pedVal);
        }

        if (mlatDeg < hallH0) {
            hallConductances[i] = std::max(0.0, hallVal);
        } else {
            hallConductances[i] = std::max(0.55, hallVal);
        }
    }
    return _auroralConductance;
}

double HardyConductance::_fourierSeries(nlohmann::json coefficients,
                                        double mlt) {
    using std::cos, std::sin;
    double coefficient = coefficients["const"].get<double>();
    for (int i = 0; i < 6; i++) {
        coefficient += coefficients["cos"][i].get<double>() *
                           cos((i + 1) * mlt * M_PI / 12) +
                       coefficients["sin"][i].get<double>() *
                           sin((i + 1) * mlt * M_PI / 12);
    }
    return coefficient;
}

double HardyConductance::_epsteinFunction(double h, double r, double h0,
                                          double S1, double S2) {
    return r + S1 * (h - h0) +
           (S2 - S1) * std::log((1.0 - (S1 / (S2 * std::exp(-(h - h0))))) /
                                (1.0 - (S1 / S2)));
}
