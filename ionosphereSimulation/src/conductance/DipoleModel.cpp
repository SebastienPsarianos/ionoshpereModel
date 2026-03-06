#include "ionosphere/conductance/DipoleModel.hpp"

#include <cmath>

DipoleModel::DipoleModel() {
    using Eigen::Vector3d;

    // TODO: Read this from table and cite
    double g_10 = -29404.8;
    double g_11 = -1450.9;
    double h_11 = 4652.5;

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

double DipoleModel::computeMLT(MagGeo subsolar, MagGeo observer) const {
    double mlt =
        (observer.longitude - subsolar.longitude) * (180.0 / M_PI) / 15.0 +
        12.0;
    return std::fmod(mlt + 24.0, 24.0);
}

const Eigen::Matrix3d& DipoleModel::rotationMatrix() const {
    return _rotationMatrix;
}
