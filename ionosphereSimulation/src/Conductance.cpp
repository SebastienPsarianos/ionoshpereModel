#include "Conductance.hpp"
#include <cmath>

void Conductance::calculateConductance(
    std::shared_ptr<GridSet<Coords>> coords,
    std::shared_ptr<GridSet<Sigma>> sigma,
    std::shared_ptr<GridSet<HppSigma>> hppSigma,
    std::shared_ptr<GridSet<DSigma>> dSigma, int nTh, int nPh) {
    _calcSigma(coords, sigma, hppSigma, nTh, nPh);
    _calcSigmaDer(coords, sigma, dSigma, nTh, nPh);
}

void Conductance::_calcSigma(std::shared_ptr<GridSet<Coords>> coords,
                             std::shared_ptr<GridSet<Sigma>> sigma,
                             std::shared_ptr<GridSet<HppSigma>> hppSigma,
                             int nTh, int nPh) {
    for (int th = 0; th < nTh; th++) {
        for (int ph = 0; ph < nPh; ph++) {
            double cos = std::cos((*coords)(Coords::TH, th, ph));
            double cos2 = std::pow(cos, 2);
            double sin = std::sin((*coords)(Coords::TH, th, ph));
            double sin2 = std::pow(sin, 2);
            double tan = sin / cos;
            double cot2 = std::pow(cos / sin, 2);

            (*sigma)(Sigma::THTH, th, ph) =
                ((*hppSigma)(HppSigma::PARALLEL, th, ph) *
                 (*hppSigma)(HppSigma::PEDERSON, th, ph) * (1 + 3 * cos2)) /
                ((*hppSigma)(HppSigma::PARALLEL, th, ph) * cos2 +
                 (*hppSigma)(HppSigma::PEDERSON, th, ph) * sin2);

            (*sigma)(Sigma::THPH, th, ph) =
                (2 * (*hppSigma)(HppSigma::PARALLEL, th, ph) *
                 (*hppSigma)(HppSigma::PEDERSON, th, ph) *
                 std::sqrt(1 + 3 * cos2)) /
                (4 * (*hppSigma)(HppSigma::PARALLEL, th, ph) * cos +
                 (*hppSigma)(HppSigma::PEDERSON, th, ph) * sin * tan);

            (*sigma)(Sigma::PHPH, th, ph) =
                (*hppSigma)(HppSigma::PEDERSON, th, ph) +
                std::pow((*hppSigma)(HppSigma::HALL, th, ph), 2) /
                    (4 * (*hppSigma)(HppSigma::PARALLEL, th, ph) * cot2 +
                     (*hppSigma)(HppSigma::PEDERSON, th, ph));
        }
    }
}

void Conductance::_calcSigmaDer(std::shared_ptr<GridSet<Coords>> coords,
                                std::shared_ptr<GridSet<Sigma>> sigma,
                                std::shared_ptr<GridSet<DSigma>> dSigma,
                                int nTh, int nPh) {
    double dTh = (*coords)(Coords::TH, 1, 0) - (*coords)(Coords::TH, 0, 0);
    double dPh = (*coords)(Coords::TH, 0, 1) - (*coords)(Coords::TH, 0, 0);

    // TODO: Theta boundary conditions

    for (int th = 1; th < nTh - 1; th++) {
        // Boundary points for phi should be continuous
        // and wrap around to the start of the phi grid.
        // So we calculate based on the adjacent left point (nPh - 2)
        // Then we set the end of the grid to be equal to the start
        (*dSigma)(DSigma::DTHPH_PH, th, 0) =
            ((*sigma)(Sigma::THPH, th, 1) -
             (*sigma)(Sigma::THPH, th, nPh - 2)) /
            (2 * dPh);
        (*dSigma)(DSigma::DTHPH_PH, th, nPh - 1) =
            (*dSigma)(DSigma::DTHPH_PH, th, 0);

        (*dSigma)(DSigma::DPHPH_PH, th, 0) =
            ((*sigma)(Sigma::PHPH, th, 1) -
             (*sigma)(Sigma::PHPH, th, nPh - 2)) /
            (2 * dPh);
        (*dSigma)(DSigma::DPHPH_PH, th, nPh - 1) =
            (*dSigma)(DSigma::DPHPH_PH, th, 0);

        // For the derivatives in the theta direction, we can calculate normally
        (*dSigma)(DSigma::DTHTH_TH, th, 0) =
            ((*sigma)(Sigma::THTH, th + 1, 0) -
             (*sigma)(Sigma::THTH, th - 1, 0)) /
            (2 * dTh);

        (*dSigma)(DSigma::DTHPH_TH, th, 0) =
            ((*sigma)(Sigma::THPH, th + 1, 0) -
             (*sigma)(Sigma::THPH, th - 1, 0)) /
            (2 * dTh);

        for (int ph = 1; ph < nPh - 1; ph++) {
            (*dSigma)(DSigma::DTHTH_TH, th, ph) =
                ((*sigma)(Sigma::THTH, th + 1, ph) -
                 (*sigma)(Sigma::THTH, th - 1, ph)) /
                (2 * dTh);
            (*dSigma)(DSigma::DTHPH_PH, th, ph) =
                ((*sigma)(Sigma::THPH, th, ph + 1) -
                 (*sigma)(Sigma::THPH, th, ph - 1)) /
                (2 * dPh);
            (*dSigma)(DSigma::DTHPH_TH, th, ph) =
                ((*sigma)(Sigma::THPH, th + 1, ph) -
                 (*sigma)(Sigma::THPH, th - 1, ph)) /
                (2 * dTh);
            (*dSigma)(DSigma::DPHPH_PH, th, ph) =
                ((*sigma)(Sigma::PHPH, th, ph + 1) -
                 (*sigma)(Sigma::PHPH, th, ph - 1)) /
                (2 * dPh);
        }
    }
}
