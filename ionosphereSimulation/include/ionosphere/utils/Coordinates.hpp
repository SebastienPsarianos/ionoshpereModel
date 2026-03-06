#pragma once
#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/utils/Grid.hpp"

class Coordinates {
  public:
    Coordinates(Ionosphere::MultiVectorRCP coords, Ionosphere::GlobalOrd nTh,
                Ionosphere::GlobalOrd nPh, Ionosphere::Scalar dTh,
                Ionosphere::Scalar dPh);

    GeoSph geoSph(Ionosphere::LocalOrd i) const;
    GeoGeo geoGeo(Ionosphere::LocalOrd i) const;

    Ionosphere::MultiVectorRCP multiVector() const;

    static GeoGeo toGeoGeo(const GeoSph& sph);
    static GeoSph toGeoSph(const GeoGeo& geo);
    static MagGeo toMagGeo(const MagSph& sph);

    const Ionosphere::GlobalOrd nTh;
    const Ionosphere::GlobalOrd nPh;
    const Ionosphere::Scalar dTh;
    const Ionosphere::Scalar dPh;

  private:
    Ionosphere::MultiVectorRCP _coords;
};

void calculateCoords(Grid<GeoSph>& coordinates);
