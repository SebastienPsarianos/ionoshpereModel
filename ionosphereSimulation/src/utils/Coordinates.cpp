#include "ionosphere/utils/Coordinates.hpp"

void Coordinates::calculateCoords(Grid<ThPh>& coordinates) {
    const size_t nTh = coordinates.nTh;
    const size_t nPh = coordinates.nPh;
    const double dPh = 2 * M_PI / nPh;
    const double dTh = M_PI / nTh;

    for (size_t th = 0; th < nTh; th++) {
        for (size_t ph = 0; ph < nPh; ph++) {
            coordinates(th, ph).th = th * dTh;
            coordinates(th, ph).ph = ph * dPh;
        }
    }
}
