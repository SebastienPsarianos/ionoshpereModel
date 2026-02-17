#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/conductance/utils.hpp"
#include "ionosphere/utils/Grid.hpp"

#include <cmath>
#include <exception>
#include <fstream>
#include <nlohmann/json>
#include <string>

EmpConductance::EmpConductance(size_t nTh, size_t nPh, double sig0, int year,
                               int day, int month, int hour)
    : _nTh(nTh), _nPh(nPh), _sig0(sig0), _hppSigma(nTh, nPh),
      _euvConductance(nTh, nPh), _auroralConductance(nTh, nPh) {
    computeGrenaTimescales(_utTime, _ttTime, year, month, day, hour);
}

int EmpConductance::computeConductance(Grid<Sigma>& sigma,
                                       Grid<DSigma>& dConductance,
                                       Grid<GeoSph>& coords) {

    // TODO: Set kp and f107 value dynamically
    if (_computeAuroralConductance(1, coords) < 0) {
        return -1;
    }
    _computeEuvConductance(175.9, coords);

    std::ofstream conductanceFile("../data/hppConductance.txt");
    if (!conductanceFile.is_open())
        throw std::exception();
    _computeHppConductance();
    _hppSigma.printWithCoords(conductanceFile, coords);
    conductanceFile.close();

    _calcSigma(sigma, coords);
    std::ofstream conductanceFileSph("../data/SphConductance.txt");
    if (!conductanceFileSph.is_open())
        throw new std::exception();

    sigma.printWithCoords(conductanceFileSph, coords);

    _calcSigmaDer(dConductance, sigma, coords);
    std::ofstream dConductanceFile("../data/dConductance.txt");
    if (!dConductanceFile.is_open())
        throw new std::exception();

    dConductance.printWithCoords(dConductanceFile, coords);
    conductanceFile.close();

    return 0;
}

void EmpConductance::_computeEuvConductance(double f107, Grid<GeoSph>& coords) {
    using std::cos, std::sin, std::sqrt;
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            // Compute Solar Zenith for this longitude / latitude using Grena
            // (2012)

            double sza =
                computeSolarZenith(_utTime, _ttTime, convert(coords(th, ph)));

            if (sza >= M_PI / 2) {
                _euvConductance(th, ph).hall = 0;
                _euvConductance(th, ph).pederson = 0;
            } else {
                _euvConductance(th, ph).hall =
                    pow(f107, 0.53) * (0.81 * cos(sza) + 0.54 * sqrt(cos(sza)));
                _euvConductance(th, ph).pederson =
                    pow(f107, 0.49) * (0.34 * cos(sza) + 0.93 * sqrt(cos(sza)));
            }
        }
    }
}

void EmpConductance::_calcSigma(Grid<Sigma>& sigma, Grid<GeoSph>& coords) {
    double cos, sin, cos2, sin2, cos3, cos4, C;
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            cos = std::cos(coords(th, ph).theta);
            sin = std::sin(coords(th, ph).theta);
            cos2 = std::pow(cos, 2);
            sin2 = std::pow(sin, 2);
            cos3 = 1.00 + 3.00 * cos2;
            cos4 = sqrt(cos3);

            C = 4.00 * _hppSigma(th, ph).parallel * cos2 +
                _hppSigma(th, ph).pederson * sin2;

            sigma(th, ph).thth = _hppSigma(th, ph).parallel *
                                 _hppSigma(th, ph).pederson * cos3 / C;
            sigma(th, ph).thph = 2.00 * _hppSigma(th, ph).parallel *
                                 _hppSigma(th, ph).hall * cos * cos4 / C;
            sigma(th, ph).phph =
                _hppSigma(th, ph).pederson +
                _hppSigma(th, ph).hall * _hppSigma(th, ph).hall * sin2 / C;
        }
    }
}

int EmpConductance::_computeAuroralConductance(int kp, Grid<GeoSph>& coords) {
    using nlohmann::json;
    if (kp > 6 || kp < 0)
        return -1;

    std::fstream pedCoeff("../src/conductance/data/hardyPedersonCoeff.json");
    if (pedCoeff) {
        json data = json::parse(pedCoeff)["K" + std::to_string(kp)];
        json fourier_max_value = data["max_value"];
        json fourier_max_latitude = data["max_latitude"];
        json fourier_up_slope = data["up_slope"];
        json fourier_down_slope = data["down_slope"];

        for (size_t th = 0; th < _nTh; th++) {
            for (size_t ph = 0; ph < _nPh; ph++) {
                MagSph observerPositionMag;
                geoCentricToDipole(observerPositionMag, coords(th, ph));

                GeoGeo subsolarPosition;
                computeSubSolar(subsolarPosition, _utTime, _ttTime);
                MagSph subsolarPositionMag;
                geoCentricToDipole(subsolarPositionMag,
                                   convert(subsolarPosition));

                double mlt = computeMagneticLocalTime(
                    convert(subsolarPositionMag), convert(observerPositionMag));

                double mlat = convert(observerPositionMag).latitude;
                double max_value = fourierSeries(fourier_max_value, mlt);
                double max_latitude = fourierSeries(fourier_max_latitude, mlt);
                double up_slope = fourierSeries(fourier_up_slope, mlt);
                double down_slope = fourierSeries(fourier_down_slope, mlt);

                _auroralConductance(th, ph).pederson = epsteinFunction(
                    mlat, max_value, max_latitude, up_slope, down_slope);

                _auroralConductance(th, ph).hall = observerPositionMag.phi;
            }
        }

    } else {
        return -1;
    }

    std::fstream hallCoeff("../src/conductance/data/hardyHallCoeff.json");
    if (hallCoeff) {
        json data = json::parse(hallCoeff)["K" + std::to_string(kp)];
        json fourier_max_value = data["max_value"];
        json fourier_max_latitude = data["max_latitude"];
        json fourier_up_slope = data["up_slope"];
        json fourier_down_slope = data["down_slope"];

        for (size_t th = 0; th < _nTh; th++) {
            for (size_t ph = 0; ph < _nPh; ph++) {

                MagSph observerPositionMag;
                geoCentricToDipole(observerPositionMag, coords(th, ph));

                GeoGeo subsolarPosition;
                computeSubSolar(subsolarPosition, _utTime, _ttTime);
                MagSph subsolarPositionMag;
                geoCentricToDipole(subsolarPositionMag,
                                   convert(subsolarPosition));

                double mlt = computeMagneticLocalTime(
                    convert(subsolarPositionMag), convert(observerPositionMag));

                double mlat = convert(observerPositionMag).latitude;
                double max_value = fourierSeries(fourier_max_value, mlt);
                double max_latitude = fourierSeries(fourier_max_latitude, mlt);
                double up_slope = fourierSeries(fourier_up_slope, mlt);
                double down_slope = fourierSeries(fourier_down_slope, mlt);

                _auroralConductance(th, ph).hall = epsteinFunction(
                    mlat, max_value, max_latitude, up_slope, down_slope);
            }
        }

    } else {
        return -1;
    }

    return 0;
}

void EmpConductance::_computeHppConductance() {
    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            _hppSigma(th, ph).hall = sqrt(_auroralConductance(th, ph).hall *
                                              _auroralConductance(th, ph).hall +
                                          _euvConductance(th, ph).hall *
                                              _euvConductance(th, ph).hall);
            _hppSigma(th, ph).pederson =
                sqrt(_auroralConductance(th, ph).pederson *
                         _auroralConductance(th, ph).pederson +
                     _euvConductance(th, ph).pederson *
                         _euvConductance(th, ph).pederson);

            _hppSigma(th, ph).parallel = 1000;
        }
    }
}

void EmpConductance::_calcSigmaDer(Grid<DSigma>& dConductance,
                                   Grid<Sigma>& sigma, Grid<GeoSph>& coords) {

    const double dTh = coords(1, 0).theta - coords(0, 0).theta;
    const double dPh = coords(0, 1).phi - coords(0, 0).phi;

    // TODO: Figure out pole behaviour for th
    for (size_t th = 1; th < _nTh - 1; th++) {
        // NOTE: Boundary points for phi should be continuous
        //      and wrap around to the start of the phi grid.
        //     Boundary points (phi = 0 and phi = _nPh - 1) are
        //      calculated as derivatives between phi = 1 and nPh-2
        dConductance(th, 0).dthph_ph =
            (sigma(th, 1).thph - sigma(th, _nPh - 2).thph) / (2 * dPh);
        dConductance(th, _nPh - 1).dthph_ph = dConductance(th, 0).dthph_ph;

        dConductance(th, 0).dphph_ph =
            (sigma(th, 1).phph - sigma(th, _nPh - 2).phph) / (2 * dPh);
        dConductance(th, _nPh - 1).dphph_ph = dConductance(th, 0).dphph_ph;

        // We also set the th derivatives to be equal accross the circular
        // boundary
        dConductance(th, 0).dthth_th =
            (sigma(th + 1, 0).thth - sigma(th - 1, 0).thth) / (2 * dTh);
        dConductance(th, _nPh - 1).dthth_th = dConductance(th, 0).dthth_th;

        dConductance(th, 0).dthph_th =
            (sigma(th + 1, 0).thph - sigma(th - 1, 0).thph) / (2 * dTh);
        dConductance(th, _nPh - 1).dthph_th = dConductance(th, 0).dthph_th;

        for (size_t ph = 1; ph < _nPh - 1; ph++) {
            dConductance(th, ph).dthth_th =
                (sigma(th + 1, ph).thth - sigma(th - 1, ph).thth) / (2 * dTh);

            dConductance(th, ph).dthph_ph =
                (sigma(th, ph + 1).thph - sigma(th, ph - 1).thph) / (2 * dPh);

            dConductance(th, ph).dthph_th =
                (sigma(th + 1, ph).thph - sigma(th - 1, ph).thph) / (2 * dTh);

            dConductance(th, ph).dphph_ph =
                (sigma(th, ph + 1).phph - sigma(th, ph - 1).phph) / (2 * dPh);
        }
    }
}
