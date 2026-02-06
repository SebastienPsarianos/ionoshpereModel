#include "ionosphere/conductance/Conductance.hpp"
#include "ionosphere/conductance/utils.hpp"
#include "ionosphere/utils/Grid.hpp"
#include <cmath>
#include <exception>
#include <fstream>

Conductance::Conductance(Grid<Sigma>& sigma, Grid<DSigma>& dConductance,
                         size_t nTh, size_t nPh, double sig0,
                         Grid<ThPh>& coords, int year, int month, int day,
                         int hour)
    : _nTh(nTh), _nPh(nPh), _sig0(sig0), _coords(coords), _sigma(sigma),
      _hppSigma(nTh, nPh), _dSigma(dConductance) {
    computeGrenaTimescales(_utTime, _ttTime, year, month, day, hour);
}

void Conductance::calculateCoefficients() {
    // TODO: Set f107 dynamically
    _euvConductance(175.9);

    std::ofstream conductanceFile("../data/EUVConductance.txt");

    if (!conductanceFile.is_open())
        throw new std::exception();
    _hppSigma.printWithCoords(conductanceFile, _coords);

    conductanceFile.close();
    //    _calcSigma();
    //    _calcSigmaDer();
}

void Conductance::_euvConductance(double f107) {
    using std::cos, std::sin, std::sqrt;
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            // Compute Solar Zenith for this longitude / latitude using Grena
            // (2012)
            double sza = computeSolarZenith(_utTime, _ttTime,
                                            M_PI / 2.0 - _coords(th, ph).th,
                                            _coords(th, ph).ph);

            _hppSigma(th, ph).parallel = sza;

            if (sza >= M_PI / 2) {
                _hppSigma(th, ph).hall = 0;
                _hppSigma(th, ph).pederson = 0;
            } else {
                _hppSigma(th, ph).hall =
                    pow(f107, 0.53) * (0.81 * cos(sza) + 0.54 * sqrt(cos(sza)));
                _hppSigma(th, ph).pederson =
                    pow(f107, 0.49) * (0.34 * cos(sza) + 0.93 * sqrt(cos(sza)));
            }
        }
    }
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

struct {
    double sinCoeff;
    double cosCoeff;
} typedef fourierCoeff_t;

double fourierSeries(std::vector<fourierCoeff_t> fourierCoefficients,
                     double mlt) {
    using std::cos, std::sin;

    double coefficient = 0;

    for (int i = 0; i < 6; i++) {
        coefficient +=
            fourierCoefficients[i].cosCoeff * cos((i + 1) * mlt / 12) +
            fourierCoefficients[i].sinCoeff * sin((i + 1) * mlt / 12);
    }

    return coefficient;
}

void Conductance::_calcSigmaDer() {
    const double dTh = _coords(1, 0).th - _coords(0, 0).th;
    const double dPh = _coords(0, 1).ph - _coords(0, 0).ph;

    // TODO: Figure out pole behaviour for th
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
