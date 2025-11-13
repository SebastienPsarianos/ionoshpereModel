#include "Solver.hpp"
#include <cmath>
#include <cstddef>

Solver::Solver(int nTh, int nPh, std::shared_ptr<GridSet<Coeff>> kappa,
               std::shared_ptr<GridSet<Ang>> coords,
               std::shared_ptr<Grid> radCurrent, Algorithm algorithm)
    : radCurrent(radCurrent), coords(coords), kappa(kappa),
      algorithm(algorithm), previousIteration(nTh, nPh, 0) {
    potential = std::make_shared<Grid>(nTh, nPh, 0);
    dTh = (*coords)(Ang::TH, 1, 0) - (*coords)(Ang::TH, 0, 0);
    dPh = (*coords)(Ang::PH, 0, 1) - (*coords)(Ang::TH, 0, 0);
    dTh2 = dTh * dTh;
    dPh2 = dPh * dPh;
}

std::shared_ptr<Grid> Solver::calculatePotential() {
    // NOTE: Boundary conditions:
    // 1. u(0,Phi) = 0
    // 2. u(PI,Phi) = 0
    // 3. u(Theta,0) = u(Theta,2*PI)
    // - 1,2 are already enforced by the inital value of the grid.
    //      Therefore, the 0 and nTh-1 indices are not updated by the solver
    // - 3 is enforced by the solver by only calculating the potential at
    //   the end of a sweep and then updating the first potential to be
    //   equal

    int itt = 0;
    if (algorithm == GAUSS_SEIDEL) {
        for (itt = 0; itt < MAX_ITERATION_NUM; itt++) {
            int sweepStyle = itt % 4;
            double residual = 0;
            switch (sweepStyle) {
            case 0:
                // Left to right (From the top row)
                for (size_t th = 1; th < nTh - 1; th++) {
                    for (size_t ph = 1; ph < nPh; ph++) {
                        (*potential)(th, ph) = _gaussSeidelFormula(th, ph);
                        residual += _calculateResidual(th, ph);
                    }
                    (*potential)(th, 0) = (*potential)(th, nPh - 1);
                }
                break;

            case 1:
                // Top to bottom (From the leftmost column)
                for (size_t ph = 1; ph < nPh; ph++) {
                    for (size_t th = 1; th < nTh - 1; th++) {
                        (*potential)(th, ph) = _gaussSeidelFormula(th, ph);
                        residual += _calculateResidual(th, ph);
                        if (ph == nPh - 1) {
                            (*potential)(th, 0) = (*potential)(th, ph);
                        }
                    }
                }
                break;

            case 2:
                // Right to left (From bottom row)
                for (size_t th = nTh - 1; th > 1; --th) {
                    for (size_t ph = nPh - 1; ph > 0; --ph) {
                        (*potential)(th, ph) = _gaussSeidelFormula(th, ph);
                        residual += _calculateResidual(th, ph);
                    }
                    (*potential)(th, nPh - 1) = (*potential)(th, 0);
                }
                break;

            case 3:
                // Bottom to top (From the rightmost column)
                for (size_t ph = nPh - 1; ph > 0; --ph) {
                    for (size_t th = nTh - 1; th > 1; --th) {
                        (*potential)(th, ph) = _gaussSeidelFormula(th, ph);
                        residual += _calculateResidual(th, ph);
                        if (ph == 0) {
                            (*potential)(th, nPh - 1) = (*potential)(th, ph);
                        }
                    }
                }
                break;
            }
            residual /= (nPh - 1) * (nTh - 2);
            if (residual <= RES_THRESHOLD) {
                std::cout << "Algorithm converged in " << itt
                          << " iterations with a residual of " << residual;
                return potential;
            }
        }
    }

    std::cerr << "Algorithm failed to converge in " << itt << " iterations.";
    return nullptr;
};

double Solver::_gaussSeidelFormula(size_t th, size_t ph) {
    // NOTE: Boundary conditions are enforced by calculating
    //      the last potential in the sweep as wrapping around to the other
    //      side of the grid. We do not calculate the first potential in the
    //      sweep. The calling function is responsible for updating the
    //      first potential to the value of the last potential and only
    //      calculating one potential.
    if (ph == 0) {
        return -(std::pow(std::sin((*coords)(Ang::TH, th, ph)), 2) *
                     RADIUS_EARTH * (*radCurrent)(th, ph) -
                 ((*kappa)(Coeff::THTH, th, ph) / dTh2) *
                     ((*potential)(th + 1, ph) + (*potential)(th - 1, ph)) -
                 ((*kappa)(Coeff::PHPH, th, ph) / dPh2) *
                     ((*potential)(th, ph + 1) + (*potential)(th, nPh - 2)) -
                 ((*kappa)(Coeff::TH, th, ph) / (2 * dTh)) *
                     ((*potential)(th + 1, ph) - (*potential)(th - 1, ph)) -
                 ((*kappa)(Coeff::PH, th, ph) / (2 * dPh)) *
                     ((*potential)(th, ph + 1) - (*potential)(th, nPh - 2))) /
               -((2 * (*kappa)(Coeff::THTH, th, ph)) / dTh2 +
                 (2 * (*kappa)(Coeff::PHPH, th, ph)) / dPh2);
    }
    if (ph == nPh - 1) {
        return -(std::pow(std::sin((*coords)(Ang::TH, th, ph)), 2) *
                     RADIUS_EARTH * (*radCurrent)(th, ph) -
                 ((*kappa)(Coeff::THTH, th, ph) / dTh2) *
                     ((*potential)(th + 1, ph) + (*potential)(th - 1, ph)) -
                 ((*kappa)(Coeff::PHPH, th, ph) / dPh2) *
                     ((*potential)(th, 1) + (*potential)(th, ph - 1)) -
                 ((*kappa)(Coeff::TH, th, ph) / (2 * dTh)) *
                     ((*potential)(th + 1, ph) - (*potential)(th - 1, ph)) -
                 ((*kappa)(Coeff::PH, th, ph) / (2 * dPh)) *
                     ((*potential)(th, 1) - (*potential)(th, ph - 1))) /
               -((2 * (*kappa)(Coeff::THTH, th, ph)) / dTh2 +
                 (2 * (*kappa)(Coeff::PHPH, th, ph)) / dPh2);
    }
    return -(std::pow(std::sin((*coords)(Ang::TH, th, ph)), 2) * RADIUS_EARTH *
                 (*radCurrent)(th, ph) -
             ((*kappa)(Coeff::THTH, th, ph) / dTh2) *
                 ((*potential)(th + 1, ph) + (*potential)(th - 1, ph)) -
             ((*kappa)(Coeff::PHPH, th, ph) / dPh2) *
                 ((*potential)(th, ph + 1) + (*potential)(th, ph - 1)) -
             ((*kappa)(Coeff::TH, th, ph) / (2 * dTh)) *
                 ((*potential)(th + 1, ph) - (*potential)(th - 1, ph)) -
             ((*kappa)(Coeff::PH, th, ph) / (2 * dPh)) *
                 ((*potential)(th, ph + 1) - (*potential)(th, ph - 1))) /
           -((2 * (*kappa)(Coeff::THTH, th, ph)) / dTh2 +
             (2 * (*kappa)(Coeff::PHPH, th, ph)) / dPh2);
}

double Solver::_calculateResidual(size_t th, size_t ph) {
    if (ph == 0) {
        return (*kappa)(Coeff::THTH, th, ph) *
                   ((*potential)(th + 1, ph) - 2 * (*potential)(th, ph) +
                    (*potential)(th - 1, ph) / dTh2) +
               (*kappa)(Coeff::PHPH, th, ph) *
                   ((*potential)(th, ph + 1) - 2 * (*potential)(th, ph) +
                    (*potential)(th, nPh - 2) / dPh2) +
               (*kappa)(Coeff::TH, th, ph) *
                   ((*potential)(th + 1, ph) -
                    (*potential)(th - 1, ph) / (2 * dTh)) +
               (*kappa)(Coeff::PH, th, ph) *
                   ((*potential)(th, ph + 1) -
                    (*potential)(th, nPh - 2) / (2 * dTh)) -
               std::pow(std::sin((*coords)(Ang::TH, th, ph)), 2) *
                   RADIUS_EARTH * (*radCurrent)(th, ph);
    }
    if (ph == nPh - 1) {
        return (*kappa)(Coeff::THTH, th, ph) *
                   ((*potential)(th + 1, ph) - 2 * (*potential)(th, ph) +
                    (*potential)(th - 1, ph) / dTh2) +
               (*kappa)(Coeff::PHPH, th, ph) *
                   ((*potential)(th, 1) - 2 * (*potential)(th, ph) +
                    (*potential)(th, ph - 1) / dPh2) +
               (*kappa)(Coeff::TH, th, ph) *
                   ((*potential)(th + 1, ph) -
                    (*potential)(th - 1, ph) / (2 * dTh)) +
               (*kappa)(Coeff::PH, th, ph) *
                   ((*potential)(th, 1) -
                    (*potential)(th, ph - 1) / (2 * dTh)) -
               std::pow(std::sin((*coords)(Ang::TH, th, ph)), 2) *
                   RADIUS_EARTH * (*radCurrent)(th, ph);
    }

    return (*kappa)(Coeff::THTH, th, ph) *
               ((*potential)(th + 1, ph) - 2 * (*potential)(th, ph) +
                (*potential)(th - 1, ph) / dTh2) +
           (*kappa)(Coeff::PHPH, th, ph) *
               ((*potential)(th, ph + 1) - 2 * (*potential)(th, ph) +
                (*potential)(th, ph - 1) / dPh2) +
           (*kappa)(Coeff::TH, th, ph) *
               ((*potential)(th + 1, ph) -
                (*potential)(th - 1, ph) / (2 * dTh)) +
           (*kappa)(Coeff::PH, th, ph) *
               ((*potential)(th, ph + 1) -
                (*potential)(th, ph - 1) / (2 * dTh)) -
           std::pow(std::sin((*coords)(Ang::TH, th, ph)), 2) * RADIUS_EARTH *
               (*radCurrent)(th, ph);
}
