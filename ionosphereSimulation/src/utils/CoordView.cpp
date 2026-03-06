#include "ionosphere/utils/CoordView.hpp"

#include <cmath>
#include <stdexcept>

CoordView::CoordView(Ionosphere::MultiVectorRCP coords,
                     const DipoleModel* dipole)
    : _coords(coords), _dipole(dipole) {}

GeoSph CoordView::geoSph(size_t i) const {
    auto thVals = _coords->getData(0);
    auto phVals = _coords->getData(1);
    return {.theta = thVals[i], .phi = phVals[i]};
}

GeoGeo CoordView::geoGeo(size_t i) const { return toGeoGeo(geoSph(i)); }

MagSph CoordView::magSph(size_t i) const {
    if (!_dipole)
        throw std::runtime_error(
            "CoordView: DipoleModel required for magnetic coordinates");
    return _dipole->geoCentricToDipole(geoSph(i));
}

MagGeo CoordView::magGeo(size_t i) const { return toMagGeo(magSph(i)); }

size_t CoordView::size() const {
    return static_cast<size_t>(_coords->getLocalLength());
}

GeoGeo CoordView::toGeoGeo(const GeoSph& sph) {
    GeoGeo geo;
    geo.latitude = M_PI_2 - sph.theta;
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI)
        geo.longitude -= 2.0 * M_PI;
    return geo;
}

GeoSph CoordView::toGeoSph(const GeoGeo& geo) {
    GeoSph sph;
    sph.theta = M_PI_2 - geo.latitude;
    sph.phi = geo.longitude;
    if (sph.phi < 0.0)
        sph.phi += 2.0 * M_PI;
    return sph;
}

MagGeo CoordView::toMagGeo(const MagSph& sph) {
    MagGeo geo;
    geo.latitude = M_PI_2 - sph.theta;
    geo.longitude = sph.phi;
    if (geo.longitude > M_PI)
        geo.longitude -= 2.0 * M_PI;
    return geo;
}
