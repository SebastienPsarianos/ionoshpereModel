#include "Conductance.hpp"
#include <cmath>

void Conductance::calculateConductance(GridSet<Coords>& coords,
                                       GridSet<Sigma>& sigma,
                                       GridSet<HppSigma>& hppSigma,
                                       GridSet<DSigma>& dSigma) {
    _calcSigma(coords, sigma, hppSigma);
    _calcSigmaDer(coords, sigma, dSigma);
}

void Conductance::_calcSigma(GridSet<Coords>& coords, GridSet<Sigma>& sigma,
                             GridSet<HppSigma>& hppSigma) {
    const size_t nTh = coords.nTh;
    const size_t nPh = coords.nPh;

    for (size_t th = 0; th < nTh; th++) {
        for (size_t ph = 0; ph < nPh; ph++) {
            double cos = std::cos(coords(Coords::TH, th, ph));
            double cos2 = std::pow(cos, 2);
            double sin = std::sin(coords(Coords::TH, th, ph));
            double sin2 = std::pow(sin, 2);
            double tan = sin / cos;
            double cot2 = std::pow(cos / sin, 2);

            sigma(Sigma::THTH, th, ph) =
                (hppSigma(HppSigma::PARALLEL, th, ph) *
                 hppSigma(HppSigma::PEDERSON, th, ph) * (1 + 3 * cos2)) /
                (hppSigma(HppSigma::PARALLEL, th, ph) * cos2 +
                 hppSigma(HppSigma::PEDERSON, th, ph) * sin2);

            sigma(Sigma::THPH, th, ph) =
                (2 * hppSigma(HppSigma::PARALLEL, th, ph) *
                 hppSigma(HppSigma::PEDERSON, th, ph) *
                 std::sqrt(1 + 3 * cos2)) /
                (4 * hppSigma(HppSigma::PARALLEL, th, ph) * cos +
                 hppSigma(HppSigma::PEDERSON, th, ph) * sin * tan);

            sigma(Sigma::PHPH, th, ph) =
                hppSigma(HppSigma::PEDERSON, th, ph) +
                std::pow(hppSigma(HppSigma::HALL, th, ph), 2) /
                    (4 * hppSigma(HppSigma::PARALLEL, th, ph) * cot2 +
                     hppSigma(HppSigma::PEDERSON, th, ph));
        }
    }
}

void Conductance::_calcSigmaDer(GridSet<Coords>& coords, GridSet<Sigma>& sigma,
                                GridSet<DSigma>& dSigma) {

    const size_t nTh = coords.nTh;
    const size_t nPh = coords.nPh;
    const double dTh = coords(Coords::TH, 1, 0) - coords(Coords::TH, 0, 0);
    const double dPh = coords(Coords::TH, 0, 1) - coords(Coords::TH, 0, 0);

    // TODO: Theta boundary conditions

    for (size_t th = 1; th < nTh - 1; th++) {
        // Boundary points for phi should be continuous
        // and wrap around to the start of the phi grid.
        // So we calculate based on the adjacent left point (nPh - 2)
        // Then we set the end of the grid to be equal to the start
        dSigma(DSigma::DTHPH_PH, th, 0) =
            (sigma(Sigma::THPH, th, 1) - sigma(Sigma::THPH, th, nPh - 2)) /
            (2 * dPh);
        dSigma(DSigma::DTHPH_PH, th, nPh - 1) = dSigma(DSigma::DTHPH_PH, th, 0);

        dSigma(DSigma::DPHPH_PH, th, 0) =
            (sigma(Sigma::PHPH, th, 1) - sigma(Sigma::PHPH, th, nPh - 2)) /
            (2 * dPh);
        dSigma(DSigma::DPHPH_PH, th, nPh - 1) = dSigma(DSigma::DPHPH_PH, th, 0);

        // For the derivatives in the theta direction, we can calculate normally
        dSigma(DSigma::DTHTH_TH, th, 0) =
            (sigma(Sigma::THTH, th + 1, 0) - sigma(Sigma::THTH, th - 1, 0)) /
            (2 * dTh);

        dSigma(DSigma::DTHPH_TH, th, 0) =
            (sigma(Sigma::THPH, th + 1, 0) - sigma(Sigma::THPH, th - 1, 0)) /
            (2 * dTh);

        for (size_t ph = 1; ph < nPh - 1; ph++) {
            dSigma(DSigma::DTHTH_TH, th, ph) =
                (sigma(Sigma::THTH, th + 1, ph) -
                 sigma(Sigma::THTH, th - 1, ph)) /
                (2 * dTh);
            dSigma(DSigma::DTHPH_PH, th, ph) =
                (sigma(Sigma::THPH, th, ph + 1) -
                 sigma(Sigma::THPH, th, ph - 1)) /
                (2 * dPh);
            dSigma(DSigma::DTHPH_TH, th, ph) =
                (sigma(Sigma::THPH, th + 1, ph) -
                 sigma(Sigma::THPH, th - 1, ph)) /
                (2 * dTh);
            dSigma(DSigma::DPHPH_PH, th, ph) =
                (sigma(Sigma::PHPH, th, ph + 1) -
                 sigma(Sigma::PHPH, th, ph - 1)) /
                (2 * dPh);
        }
    }
}
