#include "Conductance.hpp"
#include "Grid.hpp"
#include <cmath>
#include <fstream>

Conductance::Conductance(Grid<Coeff>& kappa, Grid<Sigma>& sigma, size_t nTh,
                         size_t nPh, double sig0, double sigP, double sigH,
                         Grid<ThPh>& coords)
    : _nTh(nTh), _nPh(nPh), _coords(coords), _sigma(sigma), _kappa(kappa),
      _hppSigma(nTh, nPh,
                HppSigma{.hall = sigH, .pederson = sigP, .parallel = sig0}),
      _dSigma(nTh, nPh) {}

void Conductance::calculateCoefficients() {
    _calcSigma();
    std::ofstream test("../data/conductance.txt");
    _sigma.printWithCoords(test, _coords);
    test.close();
    _calcSigmaDer();
    std::ofstream dsigma("../data/dconductance.txt");
    _dSigma.printWithCoords(dsigma, _coords);
    dsigma.close();
    _calcCoeff();
    std::ofstream kappa("../data/kappa.txt");
    _kappa.printWithCoords(kappa, _coords);
    kappa.close();
}

void Conductance::_calcSigma() {

    double cos, sin, cos2, sin2, cos3, cos4, C;
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            cos = std::cos(_coords(th, ph).th);
            sin = std::sin(_coords(th, ph).th);
            cos2 = std::pow(cos, 2);
            sin2 = std::pow(sin, 2);
            cos3 = 1.00 + 3.00 * cos2;
            cos4 = sqrt(cos3);
            C = 4.00 * _hppSigma(th, ph).parallel * cos2 +
                _hppSigma(th, ph).pederson * sin2;

            _sigma(th, ph).thth = _hppSigma(th, ph).parallel *
                                  _hppSigma(th, ph).pederson * cos3 / C;
            _sigma(th, ph).thph = 2.00 * _hppSigma(th, ph).parallel *
                                  _hppSigma(th, ph).hall * cos * cos4 / C;
            _sigma(th, ph).phph =
                _hppSigma(th, ph).pederson +
                _hppSigma(th, ph).hall * _hppSigma(th, ph).hall * sin2 / C;
        }
    }
}

void Conductance::_calcSigmaDer() {
    const double dTh = _coords(1, 0).th - _coords(0, 0).th;
    const double dPh = _coords(0, 1).ph - _coords(0, 0).ph;

    // NOTE: The original code didn't set th at 0 or at nTh - 1
    for (size_t th = 1; th < _nTh - 1; th++) {
        // NOTE: Boundary points for phi should be continuous
        //      and wrap around to the start of the phi grid.
        //     Boundary points (phi = 0 and phi = _nPh - 1) are
        //      calculated as derrivatives between phi = 1 and nPh-2
        _dSigma(th, 0).dthph_ph =
            (_sigma(th, 1).thph - _sigma(th, _nPh - 2).thph) / (2 * dPh);
        _dSigma(th, _nPh - 1).dthph_ph = _dSigma(th, 0).dthph_ph;

        _dSigma(th, 0).dphph_ph =
            (_sigma(th, 1).phph - _sigma(th, _nPh - 2).phph) / (2 * dPh);
        _dSigma(th, _nPh - 1).dphph_ph = _dSigma(th, 0).dphph_ph;

        // We also set the th derivatives to be equal accross the circular
        // boundary
        _dSigma(th, 0).dthth_th =
            (_sigma(th + 1, 0).thth - _sigma(th - 1, 0).thth) / (2 * dTh);
        _dSigma(th, _nPh - 1).dthth_th = _dSigma(th, 0).dthth_th;

        _dSigma(th, 0).dthph_th =
            (_sigma(th + 1, 0).thph - _sigma(th - 1, 0).thph) / (2 * dTh);
        _dSigma(th, _nPh - 1).dthph_th = _dSigma(th, 0).dthph_th;

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
