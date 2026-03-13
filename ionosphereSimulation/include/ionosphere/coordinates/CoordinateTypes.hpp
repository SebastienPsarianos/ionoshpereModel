#pragma once

namespace Ionosphere {
struct GeoSph {
    double theta; // [0, Pi] Colatitude (0 is north pole)
    double phi;   // [0, 2Pi]Longitude
};

struct GeoGeo {
    double latitude;  // [-PI/2, PI/2] (PI/2 is North Pole)
    double longitude; // [-PI, PI]
};

struct MagSph {
    double theta; // [0, Pi] Colatitude (0 is north pole)
    double phi;   // [0, 2Pi] Longitude
};

struct MagGeo {
    double latitude;  // [-PI/2, PI/2] (PI/2 is North Pole)
    double longitude; // [-PI, PI]
};

GeoGeo toGeoGeo(const GeoSph& sph);
GeoSph toGeoSph(const GeoGeo& geo);
MagGeo toMagGeo(const MagSph& sph);
MagSph toMagSph(const MagGeo& geo);

} // namespace Ionosphere
