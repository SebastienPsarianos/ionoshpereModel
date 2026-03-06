#pragma once
#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/conductance/DipoleModel.hpp"
#include "ionosphere/utils/Grid.hpp"

class CoordView {
  public:
    CoordView(Ionosphere::MultiVectorRCP coords,
              const DipoleModel* dipole = nullptr);

    GeoSph geoSph(size_t i) const;
    GeoGeo geoGeo(size_t i) const;
    MagSph magSph(size_t i) const; // requires DipoleModel
    MagGeo magGeo(size_t i) const; // requires DipoleModel
    size_t size() const;

    static GeoGeo toGeoGeo(const GeoSph& sph);
    static GeoSph toGeoSph(const GeoGeo& geo);
    static MagGeo toMagGeo(const MagSph& sph);

  private:
    Ionosphere::MultiVectorRCP _coords;
    const DipoleModel* _dipole; // non-owning, nullable
};
