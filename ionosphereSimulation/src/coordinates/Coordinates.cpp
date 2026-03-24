#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/coordinates/CoordinateTypes.hpp"

#include <cmath>

using namespace Ionosphere;

Coordinates::Coordinates(MultiVectorRCP coords,
                         Teuchos::RCP<DipoleModel> dipoleModel,
                         Teuchos::RCP<SolarModel> solarModel, GlobalOrd nTh,
                         GlobalOrd nPh, Scalar dTh, Scalar dPh)
    : dTh(dTh), dPh(dPh), nTh(nTh), nPh(nPh), _dipoleModel(dipoleModel),
      _solarModel(solarModel), _coords(coords) {
    _subsolarPoint = toMagGeo(dipoleModel->geoCentricToDipole(
        toGeoSph(solarModel->computeSubSolar())));
}

MagSph Coordinates::localIdx2MagSph(LocalOrd i) const {
    auto thVals = _coords->getData(0);
    auto phVals = _coords->getData(1);
    return {.theta = thVals[i], .phi = phVals[i]};
}

MagGeo Coordinates::localIdx2MagGeo(LocalOrd i) const {
    return toMagGeo(localIdx2MagSph(i));
}

GeoSph Coordinates::localIdx2GeoSph(LocalOrd i) const {
    return _dipoleModel->dipoleToGeoCentric(localIdx2MagSph(i));
}

GeoGeo Coordinates::localIdx2GeoGeo(LocalOrd i) const {
    return toGeoGeo(_dipoleModel->dipoleToGeoCentric(localIdx2MagSph(i)));
}

double Coordinates::localIdx2Sza(LocalOrd i) const {
    return _solarModel->computeZenith(localIdx2GeoGeo(i));
}

double Coordinates::localIdx2Mlt(LocalOrd i) const {
    double mlt = (localIdx2MagGeo(i).longitude - _subsolarPoint.longitude) *
                     (12 / M_PI) +
                 12.0;
    return std::fmod(mlt + 24.0, 24.0);
}

double Coordinates::getMlt(MagGeo geo) {
    double mlt =
        (geo.longitude - _subsolarPoint.longitude) * (12 / M_PI) + 12.0;
    return std::fmod(mlt + 24.0, 24.0);
}

GlobalOrd Coordinates::thPh2GlobalIdx(GlobalOrd thId, GlobalOrd phId) const {
    return thId + phId * nTh;
}

std::pair<GlobalOrd, GlobalOrd>
Coordinates::globalIdx2ThetaPhi(GlobalOrd globalIdx) const {
    return {globalIdx % nTh, globalIdx / nTh};
}

Ionosphere::MultiVectorRCP Coordinates::multiVector() const { return _coords; }
