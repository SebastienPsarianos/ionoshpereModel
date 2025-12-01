#include "EField.hpp"
#include "Constants.hpp"
#include <cmath>

std::shared_ptr<GridSet<Ang>>
calculateEField(std::shared_ptr<Grid> potential,
                std::shared_ptr<GridSet<Ang>> coords) {

    auto eField =
        std::make_shared<GridSet<Ang>>(potential->nTh, potential->nPh);
    double dTh = (*coords)(Ang::TH, 0, 0) - (*coords)(Ang::TH, 1, 0);
    double dPh = (*coords)(Ang::PH, 0, 0) - (*coords)(Ang::PH, 0, 1);

    for (size_t th = 1; th < potential->nTh - 1; th++) {
        for (size_t ph = 1; ph < potential->nPh - 1; ph++) {
            (*eField)(Ang::TH, th, ph) =
                (1 / IONOSPHERE_Radius_Earth) *
                ((*potential)(th + 1, ph) - (*potential)(th - 1, ph)) /
                (2 * dTh);
            (*eField)(Ang::PH, th, ph) =
                (1 / (IONOSPHERE_Radius_Earth *
                      std::sin((*coords)(Ang::TH, th, ph)))) *
                ((*potential)(th, ph + 1) - (*potential)(th, ph - 1)) /
                (2 * dPh);
        }
    }

    return eField;
}
