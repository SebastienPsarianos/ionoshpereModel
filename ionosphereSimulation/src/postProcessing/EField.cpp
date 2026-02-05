#include "ionosphere/postProcessing/EField.hpp"
#include "ionosphere/Constants.hpp"
#include "ionosphere/utils/Grid.hpp"
#include <cmath>

void calculateEField(Grid<ThPh>& eField, Grid<double>& potential,
                     Grid<ThPh>& coords) {
    double dTh = coords(0, 0).th - coords(1, 0).th;
    double dPh = coords(0, 0).ph - coords(0, 1).ph;

    for (size_t th = 1; th < potential.nTh - 1; th++) {
        for (size_t ph = 1; ph < potential.nPh - 1; ph++) {
            eField(th, ph).th =
                (1 / IONOSPHERE_Radius_Earth) *
                (potential(th + 1, ph) - potential(th - 1, ph)) / (2 * dTh);
            eField(th, ph).ph =
                (1 / (IONOSPHERE_Radius_Earth * std::sin(coords(th, ph).ph))) *
                (potential(th, ph + 1) - potential(th, ph - 1)) / (2 * dPh);
        }
    }
}
