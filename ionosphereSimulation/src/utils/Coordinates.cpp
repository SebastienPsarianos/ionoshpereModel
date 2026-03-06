#include "ionosphere/utils/Coordinates.hpp"
#include <cmath>

Coordinates::Coordinates(Ionosphere::MultiVectorRCP coords, size_t nTh,
                         size_t nPh, double dTh, double dPh)
    : nTh(nTh), nPh(nPh), dTh(dTh), dPh(dPh), _coords(coords) {}

GeoSph Coordinates::geoSph(size_t i) const {
    auto thVals = _coords->getData(0);
    auto phVals = _coords->getData(1);
    return {.theta = thVals[i], .phi = phVals[i]};
}

GeoGeo Coordinates::geoGeo(size_t i) const { return toGeoGeo(geoSph(i)); }

size_t Coordinates::size() const {
    return static_cast<size_t>(_coords->getLocalLength());
}

Ionosphere::MultiVectorRCP Coordinates::multiVector() const { return _coords; }

GeoGeo Coordinates::toGeoGeo(const GeoSph& sph) {
    GeoGeo geo;
    geo.latitude = M_PI_2 - sph.theta;
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI)
        geo.longitude -= 2.0 * M_PI;
    return geo;
}

GeoSph Coordinates::toGeoSph(const GeoGeo& geo) {
    GeoSph sph;
    sph.theta = M_PI_2 - geo.latitude;
    sph.phi = geo.longitude;
    if (sph.phi < 0.0)
        sph.phi += 2.0 * M_PI;
    return sph;
}

MagGeo Coordinates::toMagGeo(const MagSph& sph) {
    MagGeo geo;
    geo.latitude = M_PI_2 - sph.theta;
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI)
        geo.longitude -= 2.0 * M_PI;
    return geo;
}

void calculateCoords(Grid<GeoSph>& coordinates) {
    const size_t nTh = coordinates.nTh;
    const size_t nPh = coordinates.nPh;
    const double dPh = 2 * M_PI / nPh;
    const double dTh = M_PI / nTh;

    for (size_t th = 0; th < nTh; th++) {
        for (size_t ph = 0; ph < nPh; ph++) {
            coordinates(th, ph).theta = th * dTh;
            coordinates(th, ph).phi = ph * dPh;
        }
    }
}
