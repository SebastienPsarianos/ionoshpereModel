#include "ionosphere/coordinates/DipoleModel.hpp"
#include "nlohmann/json"

#include <cmath>
#include <fstream>

using namespace Ionosphere;

DipoleModel::DipoleModel(CommRCP comm, int year, int month, int day,
                         double hour) {
    using Eigen::Vector3d;
    using nlohmann::json;
    using std::string;
    using std::to_string;

    const int myRank = comm->getRank();
    const int rootRank = 0;
    std::string serializedData;
    int dataLength = -1;
    json igrfJson;

    // Calculate the interpolated IGRF coefficient based on date
    if (myRank == rootRank) {
        std::ifstream igrfFile("data/IGRFCoeff.json");

        if (!igrfFile.is_open()) {
            throw std::runtime_error(
                "Unable to open IGRF coefficient file. Should be located in "
                "data/IGRFCoeff.json");
        }

        igrfFile >> igrfJson;

        serializedData = igrfJson.dump();
        dataLength = static_cast<int>(serializedData.size());
    }

    Teuchos::broadcast(*comm, rootRank, 1, &dataLength);

    if (dataLength < 0)
        throw std::runtime_error("JSON parse error");

    if (myRank != rootRank) {
        serializedData.resize(dataLength);
    }

    Teuchos::broadcast(*comm, rootRank, dataLength, &serializedData[0]);
    igrfJson = json::parse(serializedData);

    double fractionalYear = static_cast<double>(year) +
                            static_cast<double>(month - 1) / 12.0 +
                            static_cast<double>(day) / 365.0 + hour / 8760.0;
    ;

    int epoch = year - year % 5;
    string epochString = to_string(epoch);
    string nextEpochString = to_string(epoch + 5);
    bool useSecularVariation = igrfJson.find(nextEpochString) == igrfJson.end();

    double distanceFromEpoch = fractionalYear - epoch;
    double g_10, g_11, h_11;

    if (useSecularVariation) {
        // Here we use the provided secular variation for the latest epoch
        g_10 = static_cast<double>(igrfJson[epochString]["g10"]) +
               static_cast<double>(igrfJson["sv"]["g10"]) * distanceFromEpoch;

        g_11 = static_cast<double>(igrfJson[epochString]["g11"]) +
               static_cast<double>(igrfJson["sv"]["g11"]) * distanceFromEpoch;

        h_11 = static_cast<double>(igrfJson[epochString]["h11"]) +
               static_cast<double>(igrfJson["sv"]["h11"]) * distanceFromEpoch;
    } else {
        // Here the interpolation is done between epochs
        g_10 = static_cast<double>(igrfJson[epochString]["g10"]) +
               distanceFromEpoch *
                   (static_cast<double>(igrfJson[nextEpochString]["g10"]) -
                    static_cast<double>(igrfJson[epochString]["g10"])) /
                   5.0;

        g_11 = static_cast<double>(igrfJson[epochString]["g11"]) +
               distanceFromEpoch *
                   (static_cast<double>(igrfJson[nextEpochString]["g11"]) -
                    static_cast<double>(igrfJson[epochString]["g11"])) /
                   5.0;

        h_11 = static_cast<double>(igrfJson[epochString]["h11"]) +
               distanceFromEpoch *
                   (static_cast<double>(igrfJson[nextEpochString]["h11"]) -
                    static_cast<double>(igrfJson[epochString]["h11"])) /
                   5.0;
    }

    Vector3d z_geo(0, 0, 1);

    // Calculate centered dipole basis
    Vector3d z_cd = -Vector3d(g_11, h_11, g_10);
    Vector3d y_cd = z_geo.cross(z_cd);
    Vector3d x_cd = y_cd.cross(z_cd);

    // Build rotation matrix
    _rotationMatrix.row(0) = x_cd.normalized();
    _rotationMatrix.row(1) = y_cd.normalized();
    _rotationMatrix.row(2) = z_cd.normalized();
}

GeoSph DipoleModel::dipoleToGeoCentric(MagSph magCoords) const {
    using Eigen::Vector3d;
    using std::acos, std::atan2, std::cos, std::sin;

    // Convert dipole spherical to dipole Cartesian (unit vector)
    Vector3d r_cd(sin(magCoords.theta) * cos(magCoords.phi),
                  sin(magCoords.theta) * sin(magCoords.phi),
                  cos(magCoords.theta));

    // R is orthogonal, so R^{-1} = R^T
    Vector3d r_cart = _rotationMatrix.transpose() * r_cd;

    GeoSph geoCoords;
    geoCoords.theta = acos(r_cart.z());
    geoCoords.phi = atan2(r_cart.y(), r_cart.x());
    if (geoCoords.phi < 0.0) {
        geoCoords.phi += 2.0 * M_PI;
    }
    return geoCoords;
}

MagSph DipoleModel::geoCentricToDipole(GeoSph geoCoords) const {
    using Eigen::Vector3d;
    using std::acos, std::atan2, std::cos, std::sin;

    Vector3d r_cart(sin(geoCoords.theta) * cos(geoCoords.phi),
                    sin(geoCoords.theta) * sin(geoCoords.phi),
                    cos(geoCoords.theta));

    Vector3d r_cd = _rotationMatrix * r_cart;

    MagSph magCoords;
    magCoords.theta = acos(r_cd.z());
    magCoords.phi = atan2(r_cd.y(), r_cd.x());
    if (magCoords.phi < 0.0) {
        magCoords.phi += 2.0 * M_PI;
    }
    return magCoords;
}
