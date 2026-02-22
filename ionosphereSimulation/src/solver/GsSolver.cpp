#include "ionosphere/solver/GsSolver.hpp"
#include "ionosphere/Constants.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>

GsSolver::GsSolver(Grid<double>& potential, size_t nTh, size_t nPh,
                   Grid<GeoSph>& coords, Grid<double>& radCurrent,
                   Grid<Sigma>& conductance, Grid<DSigma>& dConductance)
    : _nTh(nTh), _nPh(nPh), _radCurrent(radCurrent), _coords(coords),
      _potential(potential), _conductance(conductance),
      _dConductance(dConductance), _kappa(nTh, nPh) {

    // Set up the coefficient values
    _calcCoeff();

    // Initial guess
    for (size_t th = 1; th < _nTh - 1; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            _potential(th, ph) = _radCurrent(th, ph);
            if (ph == _nPh - 1) {
                _potential(th, ph) = _radCurrent(th, 0);
            }
        }
    }

    _dTh = _coords(1, 0).theta - _coords(0, 0).theta;
    _dPh = _coords(0, 1).phi - _coords(0, 0).phi;
    _dTh2 = _dTh * _dTh;
    _dPh2 = _dPh * _dPh;
}

void GsSolver::calculatePotential() {
    // NOTE: Boundary conditions:
    // 1. u(0,Phi) = 0
    // 2. u(Pi,Phi) = u(Pi-dTh, Phi) + dTh * J_r(Pi, Phi) * sin2(Pi) * R^2 /
    // SigmaThTH(Pi, Phi)
    // 3. u(Theta,0) = u(Theta,2*PI)
    // - 1,2 are already enforced by the inital value of the grid.
    //      Therefore, the 0 and _nTh-1 indices are not updated by the solver
    // - 3 is enforced by the solver by only calculating the potential at
    //   the end of a sweep and then updating the first potential to be
    //   equal

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

        // NOTE: Enforce theta boundary conditions.  Need to determine
        //          1. Is the sin even required since (sin(pi) = 0)
        //          2. Why are we using this boundary condition?
        for (size_t ph = 0; ph < _nPh; ph++) {
            double sin2 = std::sin(_coords(_nTh - 1, ph).theta) *
                          std::sin(_coords(_nTh - 1, ph).theta);

            _potential(_nTh - 1, ph) =
                _potential(_nTh - 2, ph) + _dTh * _radCurrent(_nTh - 1, ph) *
                                               sin2 * RADIUS_EARTH_2 /
                                               _conductance(_nTh - 1, ph).thth;

            // TODO: Don't do this every iteration
            _potential(0, ph) = 0;
        }

        // Red black ordering
        for (size_t th = 1; th < _nTh - 1; th++) {
            for (size_t ph = (itt + th) % 2; ph < _nPh; ph = ph + 2) {
                uOld = _potential(th, ph);
                uNew = _gaussSeidelFormula(th, ph);

                _potential(th, ph) = uNew;
                gridDiff = std::abs(uOld - uNew);
                if (gridDiff > maxDiff)
                    maxDiff = gridDiff;

                gridNorm += std::abs(uNew);
            }
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
};

double GsSolver::_gaussSeidelFormula(size_t th, size_t ph) {
    // NOTE: Boundary conditions are enforced by calculating
    //      all potentials at ph = 0 and _nPh - 1 as functions of
    //      ph = 1 and _nPh-2. This will ensure continuity in the solution.

    size_t right = ph == _nPh - 1 ? 1 : ph + 1;
    size_t left = ph == 0 ? _nPh - 2 : ph - 1;

    // TODO: Figure out if it would be useful to calculate these values only
    // once
    double sin = std::sin(_coords(th, ph).theta);
    double sin2 = sin * sin;

    return ((_kappa(th, ph).thth / _dTh2) *
                (_potential(th + 1, ph) + _potential(th - 1, ph)) +

            (_kappa(th, ph).phph / _dPh2) *
                (_potential(th, right) + _potential(th, left)) +

            (_kappa(th, ph).th / (2 * _dTh)) *
                (_potential(th + 1, ph) - _potential(th - 1, ph)) +

            (_kappa(th, ph).ph / (2 * _dPh)) *
                (_potential(th, right) - _potential(th, left)) -
            sin2 * RADIUS_EARTH_2 * _radCurrent(th, ph)) /

           ((2 * _kappa(th, ph).thth) / _dTh2 +
            (2 * _kappa(th, ph).phph) / _dPh2);
}

double GsSolver::_calculateResidual() {
    double res = 0;
    double norm = 0;

    for (size_t th = 1; th < _nTh - 1; th++) {
        for (size_t ph = 1; ph < _nPh; ph++) {
            size_t right = ph == _nPh - 1 ? 1 : ph + 1;
            double sin = std::sin(_coords(th, ph).theta);
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

void GsSolver::_calcCoeff() {
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            double cos = std::cos(_coords(th, ph).theta);
            double sin = std::sin(_coords(th, ph).theta);
            double sin2 = std::pow(sin, 2);

            _kappa(th, ph).thth = sin2 * _conductance(th, ph).thth;
            _kappa(th, ph).phph = _conductance(th, ph).phph;
            _kappa(th, ph).th = sin2 * _dConductance(th, ph).dthth_th +
                                sin * cos * _conductance(th, ph).thth -
                                sin * _dConductance(th, ph).dthph_ph;
            _kappa(th, ph).ph = sin * _dConductance(th, ph).dthph_th +
                                _dConductance(th, ph).dphph_ph;
        }
    }
}
