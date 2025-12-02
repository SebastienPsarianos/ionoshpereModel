#include "Conductance.hpp"
#include "Grid.hpp"
#include <cmath>
#include <fstream>

Conductance::Conductance(Grid<Coeff>& kappa, size_t nTh, size_t nPh,
                         double sig0, double sigP, double sigH,
                         Grid<ThPh>& coords)
    : _nTh(nTh), _nPh(nPh), _coords(coords), _sigma(nTh, nPh),
      _hppSigma(nTh, nPh, HppSigma{sig0, sigP, sigH}), _dSigma(nTh, nPh),
      _kappa(kappa) {}

void Conductance::calculateCoefficients() {
    _calcSigma();
    std::ofstream test("../data/conductance.txt");
    _sigma.printWithCoords(test, _coords);
    test.close();

    _calcSigmaDer();
    _calcCoeff();
}

void Conductance::_calcSigma() {
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            double cos = std::cos(_coords(th, ph).th);
            double cos2 = std::pow(cos, 2);
            double sin = std::sin(_coords(th, ph).th);
            double sin2 = std::pow(sin, 2);
            double tan = sin / cos;
            double cot2 = std::pow(cos / sin, 2);

            _sigma(th, ph).thth =
                (_hppSigma(th, ph).parallel * _hppSigma(th, ph).pederson *
                 (1 + 3 * cos2)) /
                (_hppSigma(th, ph).parallel * cos2 +
                 _hppSigma(th, ph).pederson * sin2);

            _sigma(th, ph).thph =
                (2 * _hppSigma(th, ph).parallel * _hppSigma(th, ph).pederson *
                 std::sqrt(1 + 3 * cos2)) /
                (4 * _hppSigma(th, ph).parallel * cos +
                 _hppSigma(th, ph).pederson * sin * tan);

            _sigma(th, ph).phph = _hppSigma(th, ph).pederson +
                                  std::pow(_hppSigma(th, ph).hall, 2) /
                                      (4 * _hppSigma(th, ph).parallel * cot2 +
                                       _hppSigma(th, ph).pederson);
        }
    }
}

void Conductance::_calcSigmaDer() {
    const double dTh = _coords(1, 0).th - _coords(0, 0).th;
    const double dPh = _coords(0, 1).ph - _coords(0, 0).ph;

    // TODO: Theta boundary conditions
    for (size_t th = 1; th < _nTh - 1; th++) {
        // Boundary points for phi should be continuous
        // and wrap around to the start of the phi grid.
        // So we calculate based on the adjacent left point (_nPh - 2)
        // Then we set the end of the grid to be equal to the start
        _dSigma(th, 0).dthph_ph =
            (_sigma(th, 1).thph - _sigma(th, _nPh - 2).thph) / (2 * dPh);
        _dSigma(th, _nPh - 1).dthph_ph = _dSigma(th, 0).dthph_ph;

        _dSigma(th, 0).dphph_ph =
            (_sigma(th, 1).phph - _sigma(th, _nPh - 2).phph) / (2 * dPh);
        _dSigma(th, _nPh - 1).dphph_ph = _dSigma(th, 0).dphph_ph;

        // For the derivatives in the theta direction, we can calculate normally
        _dSigma(th, 0).dthth_th =
            (_sigma(th + 1, 0).thth - _sigma(th - 1, 0).thth) / (2 * dTh);

        _dSigma(th, 0).dthph_th =
            (_sigma(th + 1, 0).thph - _sigma(th - 1, 0).thph) / (2 * dTh);

        for (size_t ph = 1; ph < _nPh - 1; ph++) {

            _dSigma(th, ph).dthth_th =
                (_sigma(th + 1, ph).thth - _sigma(th - 1, ph).thth) / (2 * dTh);
            _dSigma(th, ph).dthph_ph =
                (_sigma(th, ph + 1).thph - _sigma(th, ph - 1).thph) / (2 * dPh);
            _dSigma(th, ph).dthph_th =
                (_sigma(th + 1, ph).thph - _sigma(th - 1, ph).thph) / (2 * dTh);
            _dSigma(th, ph).dphph_ph =
                (_sigma(th, ph + 1).phph - _sigma(th, ph - 1).phph) / (2 * dPh);
        }
    }
}

void Conductance::_calcCoeff() {
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            double cos = std::cos(_coords(th, ph).th);
            double sin = std::sin(_coords(th, ph).th);
            double sin2 = std::pow(sin, 2);

            _kappa(th, ph).thth = sin2 * _sigma(th, ph).thth;
            _kappa(th, ph).phph = _sigma(th, ph).phph;
            _kappa(th, ph).th = sin2 * _dSigma(th, ph).dthth_th +
                                sin * cos * _sigma(th, ph).thth -
                                sin * _dSigma(th, ph).dthph_ph;
            _kappa(th, ph).ph =
                sin * _dSigma(th, ph).dthph_th + _dSigma(th, ph).dphph_ph;
        }
    }
}
