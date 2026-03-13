#include "ionosphere/coordinates/DipoleModel.hpp"

#include <cmath>

DipoleModel::DipoleModel() {
    using Eigen::Vector3d;

    // TODO:
    // - Add some way to read this from a file.
    // - Add the interpolation parameter.
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
