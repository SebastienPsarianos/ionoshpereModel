#include <cmath>
#include <ionosphere/coordinates/CoordinateTypes.hpp>

Ionosphere::GeoGeo Ionosphere::toGeoGeo(const Ionosphere::GeoSph& sph) {
    GeoGeo geo;
    geo.latitude = M_PI_2 - sph.theta;
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI)
        geo.longitude -= 2.0 * M_PI;
    return geo;
};
Ionosphere::GeoSph Ionosphere::toGeoSph(const Ionosphere::GeoGeo& geo) {
    GeoSph sph;
    sph.theta = M_PI_2 - geo.latitude;
    sph.phi = geo.longitude;
    if (sph.phi < 0.0)
        sph.phi += 2.0 * M_PI;
    return sph;
}
Ionosphere::MagGeo Ionosphere::toMagGeo(const Ionosphere::MagSph& sph) {
    MagGeo geo;
    geo.latitude = M_PI_2 - sph.theta;
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI)
        geo.longitude -= 2.0 * M_PI;
    return geo;
}
Ionosphere::MagSph Ionosphere::toMagSph(const Ionosphere::MagGeo& geo) {
    MagSph sph;
    sph.theta = M_PI_2 - geo.latitude;
    sph.phi = geo.longitude;
    if (sph.phi < 0.0)
        sph.phi += 2.0 * M_PI;
    return sph;
}
