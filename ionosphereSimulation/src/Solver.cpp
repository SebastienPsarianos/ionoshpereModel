#include "Solver.hpp"
#include "Constants.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>

Solver::Solver(Grid<double>& potential, size_t nTh, size_t nPh,
               Grid<Coeff>& kappa, Grid<ThPh>& coords, Grid<double>& radCurrent,
               Grid<Sigma>& conductance, Algorithm algorithm)
    : _nTh(nTh), _nPh(nPh), _radCurrent(radCurrent), _coords(coords),
      _kappa(kappa), _potential(potential), _conductance(conductance),
      _algorithm(algorithm) {

    // Initial guess
    for (size_t th = 1; th < _nTh - 1; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            _potential(th, ph) = _radCurrent(th, ph);
            if (ph == _nPh - 1) {
                _potential(th, ph) = _radCurrent(th, 0);
            }
        }
    }

    _dTh = _coords(1, 0).th - _coords(0, 0).th;
    _dPh = _coords(0, 1).ph - _coords(0, 0).th;
    _dTh2 = _dTh * _dTh;
    _dPh2 = _dPh * _dPh;
}

void Solver::calculatePotential() {
    // NOTE: Boundary conditions:
    // 1. u(0,Phi) = 0
    // 2. u(Pi,Phi) = u(Pi-dTh, Phi) + dTh * J_r(Pi, Phi) / SigmaThTH(Pi, Phi) ?
    // 3. u(Theta,0) = u(Theta,2*PI)
    // - 1,2 are already enforced by the inital value of the grid.
    //      Therefore, the 0 and _nTh-1 indices are not updated by the solver
    // - 3 is enforced by the solver by only calculating the potential at
    //   the end of a sweep and then updating the first potential to be
    //   equal

    if (_algorithm == GAUSS_SEIDEL) {
        int itt = 0;
        double uNew;
        double uOld;
        double gridDiff;
        double gridNorm;
        double maxDiff;

        for (itt = 0; itt < MAX_ITERATION_NUM; itt++) {
            gridDiff = 0;
            maxDiff = 0;
            gridNorm = 0;
            // Red black ordering
            for (size_t th = 1; th < _nTh - 1; th++) {
                for (size_t ph = 1 + (itt + th) % 2; ph < _nPh; ph = ph + 2) {
                    uOld = _potential(th, ph);
                    uNew = _gaussSeidelFormula(th, ph);

                    _potential(th, ph) = uNew;
                    gridDiff = std::abs(uOld - uNew);
                    if (gridDiff > maxDiff)
                        maxDiff = gridDiff;

                    gridNorm += std::abs(uNew);
                }
                _potential(th, 0) = _potential(th, _nPh - 1);
            }

            // NOTE: Enforce theta boundary condition. (other pole is enforced
            // by initial grid value)
            for (size_t ph = 0; ph < _nPh; ph++) {
                _potential(_nTh - 1, ph) = _potential(_nTh - 2, ph) +
                                           _dTh * _radCurrent(_nTh - 1, ph) /
                                               _conductance(_nTh - 1, ph).thth;
            }

            gridNorm = gridNorm / (_nPh * _nTh);

            if (maxDiff / gridNorm <= RES_THRESHOLD) {
                std::cout << "Algorithm converged in " << itt
                          << " iterations with a residual of "
                          << this->_calculateResidual() << std::endl;
                return;
            }
        }
        std::cerr << "Algorithm failed to converge in " << itt
                  << " iterations, with residual " << this->_calculateResidual()
                  << std::endl;
    }
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
    // TODO: Implement system to calculate unchanging values only once
    size_t right = ph == _nPh - 1 ? 1 : ph + 1;
    double sin = std::sin(_coords(th, ph).th);
    double sin2 = sin * sin;

    return ((_kappa(th, ph).thth / _dTh2) *
                (_potential(th + 1, ph) + _potential(th - 1, ph)) +

            (_kappa(th, ph).phph / _dPh2) *
                (_potential(th, right) + _potential(th, ph - 1)) +

            (_kappa(th, ph).th / (2 * _dTh)) *
                (_potential(th + 1, ph) - _potential(th - 1, ph)) +

            (_kappa(th, ph).ph / (2 * _dPh)) *
                (_potential(th, right) - _potential(th, ph - 1)) -
            sin2 * RADIUS_EARTH_2 * _radCurrent(th, ph)) /

           ((2 * _kappa(th, ph).thth) / _dTh2 +
            (2 * _kappa(th, ph).phph) / _dPh2);
}

double Solver::_calculateResidual() {
    double res = 0;
    double norm = 0;

    for (size_t th = 1; th < _nTh - 1; th++) {
        for (size_t ph = 1; ph < _nPh; ph++) {
            size_t right = ph == _nPh - 1 ? 1 : ph + 1;
            double sin = std::sin(_coords(th, ph).th);
            double sin2 = sin * sin;
            norm += std::abs(sin2 * RADIUS_EARTH_2 * _radCurrent(th, ph));
            res +=
                std::abs((_kappa(th, ph).thth / _dTh2) *
                             (_potential(th + 1, ph) - 2 * _potential(th, ph) +
                              _potential(th - 1, ph)) +
                         (_kappa(th, ph).phph / _dPh2) *
                             (_potential(th, right) - 2 * _potential(th, ph) +
                              _potential(th, ph - 1)) +
                         (_kappa(th, ph).th / (2 * _dTh)) *
                             (_potential(th + 1, ph) - _potential(th - 1, ph)) +
                         (_kappa(th, ph).ph / (2 * _dPh)) *
                             (_potential(th, right) - _potential(th, ph - 1)) -
                         sin2 * RADIUS_EARTH_2 * _radCurrent(th, ph));
        }
    }
    if (norm == 0)
        return res;
    return res / norm;
}
