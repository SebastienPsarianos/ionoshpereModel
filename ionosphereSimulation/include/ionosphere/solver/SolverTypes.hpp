#pragma once

#include <string>

namespace Ionosphere {

using BCOption = enum {
    GHOSTRING,
    POLEDIVTHM,
};

using EqnOption = enum { OAMM };

inline EqnOption equationTypeFromString(const std::string& s) {
    if (s == "OAmm")
        return OAMM;
    throw std::invalid_argument("Unknown equation type: " + s);
}

inline BCOption bcTypeFromString(const std::string& s) {
    if (s == "GhostRing")
        return GHOSTRING;
    else if (s == "PoleDivThm")
        return POLEDIVTHM;
    throw std::invalid_argument("Unknown boundary condition type: " + s);
}

}; // namespace Ionosphere
