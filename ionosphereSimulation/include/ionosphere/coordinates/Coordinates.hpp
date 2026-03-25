#pragma once
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/CoordinateTypes.hpp"
#include "ionosphere/coordinates/DipoleModel.hpp"
#include "ionosphere/coordinates/SolarModel.hpp"

class Coordinates {
  public:
    Coordinates(Ionosphere::MultiVectorRCP coords,
                Teuchos::RCP<DipoleModel> dipoleModel,
                Teuchos::RCP<SolarModel> solarModel, Ionosphere::GlobalOrd nTh,
                Ionosphere::GlobalOrd nPh, Ionosphere::Scalar dTh,
                Ionosphere::Scalar dPh);

    // Coordinate Lookup
    Ionosphere::MagSph localIdx2MagSph(Ionosphere::LocalOrd i) const;
    Ionosphere::MagGeo localIdx2MagGeo(Ionosphere::LocalOrd i) const;
    Ionosphere::GeoSph localIdx2GeoSph(Ionosphere::LocalOrd i) const;
    Ionosphere::GeoGeo localIdx2GeoGeo(Ionosphere::LocalOrd i) const;
    double localIdx2Sza(Ionosphere::LocalOrd i) const;
    double localIdx2Mlt(Ionosphere::LocalOrd i) const;

    // Utilities for converting from theta and phi indexes back and forth
    // between vector indices
    Ionosphere::GlobalOrd thPh2GlobalIdx(Ionosphere::GlobalOrd thId,
                                         Ionosphere::GlobalOrd phId) const;
    std::pair<Ionosphere::GlobalOrd, Ionosphere::GlobalOrd>
    globalIdx2ThetaPhi(Ionosphere::GlobalOrd globalIdx) const;

    // Grid neighbor with the phi wrap-around.
    struct Neighbors {
        Ionosphere::GlobalOrd up, down, left, right;
    };
    Neighbors getNeighbors(Ionosphere::GlobalOrd globalIdx) const;
    bool isPole(Ionosphere::GlobalOrd globalIdx) const;

    // Exposes underlying coordinates
    Ionosphere::MultiVectorRCP multiVector() const;

    double getMlt(Ionosphere::MagGeo geo);

    const Ionosphere::Scalar dTh;
    const Ionosphere::Scalar dPh;
    const Ionosphere::GlobalOrd nTh;
    const Ionosphere::GlobalOrd nPh;

  private:
    Ionosphere::MagGeo _subsolarPoint;

    Teuchos::RCP<DipoleModel> _dipoleModel;
    Teuchos::RCP<SolarModel> _solarModel;
    Ionosphere::MultiVectorRCP _coords;
};
