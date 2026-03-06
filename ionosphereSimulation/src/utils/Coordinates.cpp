#include "ionosphere/utils/Coordinates.hpp"
#include <cmath>

using namespace Ionosphere;

Coordinates::Coordinates(MultiVectorRCP coords, GlobalOrd nTh, GlobalOrd nPh,
                         Scalar dTh, Scalar dPh)
    : nTh(nTh), nPh(nPh), dTh(dTh), dPh(dPh), _coords(coords) {}

GeoSph Coordinates::geoSph(LocalOrd i) const {
    auto thVals = _coords->getData(0);
    auto phVals = _coords->getData(1);
    return {.theta = thVals[i], .phi = phVals[i]};
}

GeoGeo Coordinates::geoGeo(LocalOrd i) const { return toGeoGeo(geoSph(i)); }

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
    const GlobalOrd nTh = coordinates.nTh;
    const GlobalOrd nPh = coordinates.nPh;
    const Scalar dPh = 2 * M_PI / nPh;
    const Scalar dTh = M_PI / nTh;

    for (GlobalOrd th = 0; th < nTh; th++) {
        for (GlobalOrd ph = 0; ph < nPh; ph++) {
            coordinates(th, ph).theta = th * dTh;
            coordinates(th, ph).phi = ph * dPh;
        }
    }
}
