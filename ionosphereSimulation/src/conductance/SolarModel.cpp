#include "ionosphere/conductance/SolarModel.hpp"

#include <cmath>

SolarModel::SolarModel(int year, int month, int day, int hour) {
    // "Correcting" the month (Grena paper)
    if (month <= 2) {
        month += 12;
        year--;
    }

    // Equation 2 (Grena): Calculate UT time in days
    _utTime = floor(365.25 * (year - 2000)) + floor(30.6001 * (month + 1)) -
              floor(0.01 * year) + day + hour / 24.0 - 21958;

    // Equation 1 (Grena): Linear approximation for deltaT
    double deltaT = 96.4 + 0.00158 * _utTime;

    // Equation 3 (Grena): Convert to TT time
    _ttTime = _utTime + 1.1574e-5 * deltaT;
}

void SolarModel::_computeSunPosition(double& alpha, double& declination) const {
    using std::asin, std::atan2, std::cos, std::sin;

    double elipLongitude = -1.388803 + 1.720279216e-2 * _ttTime +
                           3.3366e-2 * sin(SUNFREQ * _ttTime - 0.06172) +
                           3.53e-4 * sin(2 * SUNFREQ * _ttTime - 0.1163);
    double earthInclination = 4.089567e-1 - 6.19e-9 * _ttTime;

    alpha =
        atan2(sin(elipLongitude) * cos(earthInclination), cos(elipLongitude));
    declination = asin(sin(elipLongitude) * sin(earthInclination));
}

double SolarModel::computeZenith(GeoGeo coords) const {
    using std::asin, std::cos, std::sin;

    double alpha, declination;
    _computeSunPosition(alpha, declination);

    double hourAngle =
        1.7528311 + 6.300388099 * _utTime + coords.longitude - alpha;
    hourAngle = std::fmod(hourAngle + M_PI, 2 * M_PI) - M_PI;

    double elevationAngle =
        asin(sin(coords.latitude) * sin(declination) +
             cos(coords.latitude) * cos(declination) * cos(hourAngle));

    return M_PI / 2 - elevationAngle;
}

GeoGeo SolarModel::computeSubSolar() const {
    double alpha, declination;
    _computeSunPosition(alpha, declination);

    GeoGeo subSolar;
    subSolar.longitude =
        std::fmod(alpha - 1.7528311 - 6.300388099 * _utTime, 2.0 * M_PI);

    if (subSolar.longitude > M_PI) {
        subSolar.longitude -= 2.0 * M_PI;
    } else if (subSolar.longitude < -M_PI) {
        subSolar.longitude += 2.0 * M_PI;
    }
    subSolar.latitude = declination;
    return subSolar;
}
