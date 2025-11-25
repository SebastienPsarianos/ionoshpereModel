#include "Solver.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>

Solver::Solver(int nTh, int nPh, std::shared_ptr<GridSet<Coeff>> kappa,
               std::shared_ptr<GridSet<Ang>> coords,
               std::shared_ptr<Grid> radCurrent, Algorithm algorithm)
    : nTh(nTh), nPh(nPh), radCurrent(radCurrent), coords(coords), kappa(kappa),
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
    double residual = 0;

    if (algorithm == GAUSS_SEIDEL) {
        for (itt = 0; itt < MAX_ITERATION_NUM; itt++) {
            residual = 0;
            for (size_t th = 1; th < nTh - 1; th++) {
                for (size_t ph = 1 + (itt + th) % 2; ph < nPh; ph = ph + 2) {
                    (*potential)(th, ph) = _gaussSeidelFormula(th, ph);
                    residual += _calculateResidual(th, ph);
                    std::cout << _calculateResidual(th, ph) << "\n";
                }
                (*potential)(th, 0) = (*potential)(th, nPh - 1);
            }

            residual /= (nPh - 1) * (nTh - 2);

            if (residual <= RES_THRESHOLD) {
                std::cout << "Algorithm converged in " << itt
                          << " iterations with a residual of " << residual
                          << std::endl;
                return potential;
            }
        }
    }

    std::cerr << "Algorithm failed to converge in " << itt
              << " iterations, with residual " << residual << std::endl;
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
        throw new std::out_of_range(
            "potential should not be calculated directly at ph = 0");
    }
    size_t right = ph == nPh - 1 ? 1 : ph + 1;
    double sin = std::sin((*coords)(Ang::TH, th, ph));
    double sin2 = sin * sin;

    return (((*kappa)(Coeff::THTH, th, ph) / dTh2) *
                ((*potential)(th + 1, ph) + (*potential)(th - 1, ph)) +

            ((*kappa)(Coeff::PHPH, th, ph) / dPh2) *
                ((*potential)(th, right) + (*potential)(th, ph - 1)) +

            ((*kappa)(Coeff::TH, th, ph) / (2 * dTh)) *
                ((*potential)(th + 1, ph) - (*potential)(th - 1, ph)) +

            ((*kappa)(Coeff::PH, th, ph) / (2 * dPh)) *
                ((*potential)(th, right) - (*potential)(th, ph - 1)) -
            sin2 * RADIUS_EARTH * (*radCurrent)(th, ph)) /

           ((2 * (*kappa)(Coeff::THTH, th, ph)) / dTh2 +
            (2 * (*kappa)(Coeff::PHPH, th, ph)) / dPh2);
}

double Solver::_calculateResidual(size_t th, size_t ph) {

    if (ph == 0) {
        throw new std::out_of_range(
            "Should not be calculating the residual at ph = 0");
    }

    size_t right = ph == nPh - 1 ? 1 : ph + 1;
    double sin = std::sin((*coords)(Ang::TH, th, ph));
    double sin2 = sin * sin;

    return std::abs(((*kappa)(Coeff::THTH, th, ph) / dTh2) *
                        ((*potential)(th + 1, ph) - 2 * (*potential)(th, ph) +
                         (*potential)(th - 1, ph)) +
                    ((*kappa)(Coeff::PHPH, th, ph) / dPh2) *
                        ((*potential)(th, right) - 2 * (*potential)(th, ph) +
                         (*potential)(th, ph - 1)) +
                    ((*kappa)(Coeff::TH, th, ph) / (2 * dTh)) *
                        ((*potential)(th + 1, ph) - (*potential)(th - 1, ph)) +
                    ((*kappa)(Coeff::PH, th, ph) / (2 * dPh)) *
                        ((*potential)(th, right) - (*potential)(th, ph - 1)) -
                    sin2 * RADIUS_EARTH * (*radCurrent)(th, ph));
}
