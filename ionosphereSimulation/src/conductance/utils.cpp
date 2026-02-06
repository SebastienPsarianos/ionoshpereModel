#include "ionosphere/conductance/utils.hpp"
#include <cmath>

double computeSolarZenith(double utTime, double ttTime, double latitude,
                          double longitude) {
    using std::atan2, std::sin, std::cos, std::asin, std::acos, std::floor;

    double elipLongitude = -1.388803 + 1.720279216e-2 * ttTime +
                           3.3366e-2 * sin(SUNFREQ * ttTime - 0.06172) +
                           3.53e-4 * sin(2 * SUNFREQ * ttTime - 0.1163);
    double earthInclination = 4.089567e-1 - 6.19e-9 * ttTime;

    double alpha =
        atan2(sin(elipLongitude) * cos(earthInclination), cos(elipLongitude));
    double declination = asin(sin(elipLongitude) * sin(earthInclination));

    double hourAngle = 1.7528311 + 6.300388099 * utTime + longitude - alpha;
    hourAngle = std::fmod(hourAngle + M_PI, 2 * M_PI) - M_PI;

    double elevationAngle =
        asin(sin(latitude) * sin(declination) +
             cos(latitude) * cos(declination) * cos(hourAngle));

    elevationAngle = elevationAngle - 4.26e-5 * cos(elevationAngle);

    return M_PI / 2 - elevationAngle;
}

void computeGrenaTimescales(double& utTime, double& ttTime, int year, int month,
                            int day, int hour) {
    // "Correcting" the month (grena paper)
    if (month <= 2) {
        month += 12;
        year--;
    }

    // Equation 2 (Grena): Calculate UT time in days
    utTime = floor(365.25 * (year - 2000)) + floor(30.6001 * (month + 1)) -
             floor(0.01 * year) + day + hour / 24.0 - 21958;

    // Equation 1 (Grena): Linear approximation for deltaT
    double deltaT = 96.4 + 0.00158 * utTime;

    // Equation 3 (Grena): Convert to TT time
    ttTime = utTime + 1.1574e-5 * deltaT;
}
