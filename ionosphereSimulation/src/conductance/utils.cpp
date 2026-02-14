#include "ionosphere/conductance/utils.hpp"

#include <Eigen/Dense>
#include <cmath>

/* BEGIN CONVERSIONS TODO: CLEANUP THIS STUFF */
GeoGeo convert(const GeoSph& sph) {
    GeoGeo geo;
    geo.latitude = M_PI_2 - sph.theta; // Colatitude to Latitude

    // Normalize Phi (0 to 2PI) into Longitude (-PI to PI)
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI) {
        geo.longitude -= 2.0 * M_PI;
    }
    return geo;
}

MagGeo convert(const MagSph& sph) {
    MagGeo geo;
    geo.latitude = M_PI_2 - sph.theta; // Colatitude to Latitude

    // Normalize Phi (0 to 2PI) into Longitude (-PI to PI)
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI) {
        geo.longitude -= 2.0 * M_PI;
    }
    return geo;
}

GeoSph convert(const GeoGeo& geo) {
    GeoSph sph;
    sph.theta = M_PI_2 - geo.latitude; // Latitude to Colatitude

    // Normalize Longitude (-PI to PI) into Phi (0 to 2PI)
    sph.phi = geo.longitude;
    if (sph.phi < 0.0) {
        sph.phi += 2.0 * M_PI;
    }
    return sph;
}

MagSph convert(const MagGeo& geo) {
    MagSph sph;
    sph.theta = M_PI_2 - geo.latitude; // Latitude to Colatitude

    // Normalize Longitude (-PI to PI) into Phi (0 to 2PI)
    sph.phi = geo.longitude;
    if (sph.phi < 0.0) {
        sph.phi += 2.0 * M_PI;
    }
    return sph;
}
/* END CONVERSIONS TODO: CLEANUP THIS STUFF */

void computeSunPosition(double& alpha, double& declination, double ttTime) {

    using std::atan2, std::sin, std::cos, std::asin;

    double elipLongitude = -1.388803 + 1.720279216e-2 * ttTime +
                           3.3366e-2 * sin(SUNFREQ * ttTime - 0.06172) +
                           3.53e-4 * sin(2 * SUNFREQ * ttTime - 0.1163);
    double earthInclination = 4.089567e-1 - 6.19e-9 * ttTime;

    alpha =
        atan2(sin(elipLongitude) * cos(earthInclination), cos(elipLongitude));
    declination = asin(sin(elipLongitude) * sin(earthInclination));
}

double computeSolarZenith(double utTime, double ttTime, GeoGeo coords) {
    using std::sin, std::cos, std::asin;

    double alpha, declination;
    computeSunPosition(alpha, declination, ttTime);

    double hourAngle =
        1.7528311 + 6.300388099 * utTime + coords.longitude - alpha;
    hourAngle = std::fmod(hourAngle + M_PI, 2 * M_PI) - M_PI;

    double elevationAngle =
        asin(sin(coords.latitude) * sin(declination) +
             cos(coords.latitude) * cos(declination) * cos(hourAngle));

    // elevationAngle = elevationAngle - 4.26e-5 * cos(elevationAngle);

    return M_PI / 2 - elevationAngle;
}

void computeSubSolar(GeoGeo& subSolar, double utTime, double ttTime) {
    double alpha, declination;
    computeSunPosition(alpha, declination, ttTime);

    subSolar.longitude =
        std::fmod(alpha - 1.7528311 - 6.300388099 * utTime, 2.0 * M_PI);

    if (subSolar.longitude > M_PI) {
        subSolar.longitude -= 2.0 * M_PI;
    } else if (subSolar.longitude < -M_PI) {
        subSolar.longitude += 2.0 * M_PI;
    }
    subSolar.latitude = declination;
}

double computeMagneticLocalTime(MagGeo subsolar, MagGeo observer) {
    double mlt =
        (observer.longitude - subsolar.longitude) * (180.0 / M_PI) / 15.0 +
        12.0;
    return std::fmod(mlt + 24.0, 24.0);
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

// TODO:
// 1. Store the matrix and just apply it when needed
// 2. Grab the values from the table
// 3. Look into interpolation
void geoCentricToDipole(MagSph& magCoords, GeoSph geoCoords) {
    using Eigen::Vector3d, Eigen::Matrix3d, std::sin, std::cos, std::acos,
        std::atan2;

    // TODO: Read this from table and cite
    double g_10 = -29404.8;
    double g_11 = -1450.9;
    double h_11 = 4652.5;

    Vector3d r_cart(sin(geoCoords.theta) * cos(geoCoords.phi),
                    sin(geoCoords.theta) * sin(geoCoords.phi),
                    cos(geoCoords.theta));
    Vector3d z_geo(0, 0, 1);

    // Calculate centered dipole basis
    Vector3d z_cd = -Vector3d(g_11, h_11, g_10);
    Vector3d y_cd = z_geo.cross(z_cd);
    Vector3d x_cd = y_cd.cross(z_cd);

    // Build rotation matrix
    Matrix3d globalToMagRotation;
    globalToMagRotation.row(0) = x_cd.normalized();
    globalToMagRotation.row(1) = y_cd.normalized();
    globalToMagRotation.row(2) = z_cd.normalized();

    Vector3d r_cd = globalToMagRotation * r_cart;

    magCoords.theta = acos(r_cd.z());
    magCoords.phi = atan2(r_cd.y(), r_cd.x());
    if (magCoords.phi < 0.0) {
        magCoords.phi += 2.0 * M_PI;
    }
}

double fourierSeries(nlohmann::json coefficients, double mlt) {
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

double epsteinFunction(double h, double r, double h0, double S1, double S2) {
    return r + S1 * (h - h0) +
           (S2 - S1) * std::log((1.0 - (S1 / (S2 * std::exp(-(h - h0))))) /
                                (1.0 - (S1 / S2)));
}
