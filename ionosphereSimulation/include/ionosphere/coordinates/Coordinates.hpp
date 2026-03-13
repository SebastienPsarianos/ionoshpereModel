#pragma once
#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/coordinates/DipoleModel.hpp"
#include "ionosphere/coordinates/Grid.hpp"
#include "ionosphere/coordinates/SolarModel.hpp"

class Coordinates {
  public:
    Coordinates(Ionosphere::MultiVectorRCP coords,
                Teuchos::RCP<DipoleModel> dipoleModel,
                Teuchos::RCP<SolarModel> solarModel, Ionosphere::GlobalOrd nTh,
                Ionosphere::GlobalOrd nPh, Ionosphere::Scalar dTh,
                Ionosphere::Scalar dPh);

    // Coordinate Lookup
    MagSph localIdx2MagSph(Ionosphere::LocalOrd i) const;
    MagGeo localIdx2MagGeo(Ionosphere::LocalOrd i) const;
    GeoSph localIdx2GeoSph(Ionosphere::LocalOrd i) const;
    GeoGeo localIdx2GeoGeo(Ionosphere::LocalOrd i) const;
    double localIdx2Sza(Ionosphere::LocalOrd i) const;
    double localIdx2Mlt(Ionosphere::LocalOrd i) const;

    // Utilities for converting from theta and phi indexes back and forth
    // between vector indices
    Ionosphere::GlobalOrd thPh2GlobalIdx(Ionosphere::GlobalOrd thId,
                                         Ionosphere::GlobalOrd phId) const;
    std::pair<Ionosphere::GlobalOrd, Ionosphere::GlobalOrd>
    globalIdx2ThetaPhi(Ionosphere::GlobalOrd globalIdx) const;

    // Utils for transforming coordinates
    static GeoGeo toGeoGeo(const GeoSph& sph);
    static GeoSph toGeoSph(const GeoGeo& geo);
    static MagGeo toMagGeo(const MagSph& sph);
    static MagSph toMagSph(const MagGeo& geo);

    // Exposes underlying coordinates
    Ionosphere::MultiVectorRCP multiVector() const;

    const Ionosphere::Scalar dTh;
    const Ionosphere::Scalar dPh;
    const Ionosphere::GlobalOrd nTh;
    const Ionosphere::GlobalOrd nPh;

  private:
    MagGeo _subsolarPoint;

    Teuchos::RCP<DipoleModel> _dipoleModel;
    Teuchos::RCP<SolarModel> _solarModel;
    Ionosphere::MultiVectorRCP _coords;
};
