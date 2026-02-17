#include "ionosphere/postProcessing/EField.hpp"
#include "ionosphere/Constants.hpp"
#include "ionosphere/utils/Grid.hpp"
#include <cmath>

void calculateEField(Grid<GeoSph>& eField, Grid<double>& potential,
                     Grid<GeoSph>& coords) {
    double dTh = coords(0, 0).theta - coords(1, 0).theta;
    double dPh = coords(0, 0).phi - coords(0, 1).phi;

    for (size_t th = 1; th < potential.nTh - 1; th++) {
        for (size_t ph = 1; ph < potential.nPh - 1; ph++) {
            eField(th, ph).theta =
                (1 / IONOSPHERE_Radius_Earth) *
                (potential(th + 1, ph) - potential(th - 1, ph)) / (2 * dTh);
            eField(th, ph).phi =
                (1 / (IONOSPHERE_Radius_Earth * std::sin(coords(th, ph).phi))) *
                (potential(th, ph + 1) - potential(th, ph - 1)) / (2 * dPh);
        }
    }
}
