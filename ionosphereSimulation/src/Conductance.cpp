#include "Conductance.hpp"
#include <cmath>

std::shared_ptr<GridSet<Coeff>> Conductance::calculateCoefficients() {
    _calcSigma();
    _calcSigmaDer();
    _calcCoeff();
    return _coefficients;
}

void Conductance::_calcSigma() {

    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            double cos = std::cos((*_coords)(Ang::TH, th, ph));
            double cos2 = std::pow(cos, 2);
            double sin = std::sin((*_coords)(Ang::TH, th, ph));
            double sin2 = std::pow(sin, 2);
            double tan = sin / cos;
            double cot2 = std::pow(cos / sin, 2);

            _sigma(Sigma::THTH, th, ph) =
                (_hppSigma(HppSigma::PARALLEL, th, ph) *
                 _hppSigma(HppSigma::PEDERSON, th, ph) * (1 + 3 * cos2)) /
                (_hppSigma(HppSigma::PARALLEL, th, ph) * cos2 +
                 _hppSigma(HppSigma::PEDERSON, th, ph) * sin2);

            _sigma(Sigma::THPH, th, ph) =
                (2 * _hppSigma(HppSigma::PARALLEL, th, ph) *
                 _hppSigma(HppSigma::PEDERSON, th, ph) *
                 std::sqrt(1 + 3 * cos2)) /
                (4 * _hppSigma(HppSigma::PARALLEL, th, ph) * cos +
                 _hppSigma(HppSigma::PEDERSON, th, ph) * sin * tan);

            _sigma(Sigma::PHPH, th, ph) =
                _hppSigma(HppSigma::PEDERSON, th, ph) +
                std::pow(_hppSigma(HppSigma::HALL, th, ph), 2) /
                    (4 * _hppSigma(HppSigma::PARALLEL, th, ph) * cot2 +
                     _hppSigma(HppSigma::PEDERSON, th, ph));
        }
    }
}

void Conductance::_calcSigmaDer() {
    const double _dTh = (*_coords)(Ang::TH, 1, 0) - (*_coords)(Ang::TH, 0, 0);
    const double dPh = (*_coords)(Ang::TH, 0, 1) - (*_coords)(Ang::TH, 0, 0);

    // TODO: Theta boundary conditions

    for (size_t th = 1; th < _nTh - 1; th++) {
        // Boundary points for phi should be continuous
        // and wrap around to the start of the phi grid.
        // So we calculate based on the adjacent left point (_nPh - 2)
        // Then we set the end of the grid to be equal to the start
        _dSigma(DSigma::DTHPH_PH, th, 0) =
            (_sigma(Sigma::THPH, th, 1) - _sigma(Sigma::THPH, th, _nPh - 2)) /
            (2 * dPh);
        _dSigma(DSigma::DTHPH_PH, th, _nPh - 1) =
            _dSigma(DSigma::DTHPH_PH, th, 0);

        _dSigma(DSigma::DPHPH_PH, th, 0) =
            (_sigma(Sigma::PHPH, th, 1) - _sigma(Sigma::PHPH, th, _nPh - 2)) /
            (2 * dPh);
        _dSigma(DSigma::DPHPH_PH, th, _nPh - 1) =
            _dSigma(DSigma::DPHPH_PH, th, 0);

        // For the derivatives in the theta direction, we can calculate normally
        _dSigma(DSigma::DTHTH_TH, th, 0) =
            (_sigma(Sigma::THTH, th + 1, 0) - _sigma(Sigma::THTH, th - 1, 0)) /
            (2 * _dTh);

        _dSigma(DSigma::DTHPH_TH, th, 0) =
            (_sigma(Sigma::THPH, th + 1, 0) - _sigma(Sigma::THPH, th - 1, 0)) /
            (2 * _dTh);

        for (size_t ph = 1; ph < _nPh - 1; ph++) {
            _dSigma(DSigma::DTHTH_TH, th, ph) =
                (_sigma(Sigma::THTH, th + 1, ph) -
                 _sigma(Sigma::THTH, th - 1, ph)) /
                (2 * _dTh);
            _dSigma(DSigma::DTHPH_PH, th, ph) =
                (_sigma(Sigma::THPH, th, ph + 1) -
                 _sigma(Sigma::THPH, th, ph - 1)) /
                (2 * dPh);
            _dSigma(DSigma::DTHPH_TH, th, ph) =
                (_sigma(Sigma::THPH, th + 1, ph) -
                 _sigma(Sigma::THPH, th - 1, ph)) /
                (2 * _dTh);
            _dSigma(DSigma::DPHPH_PH, th, ph) =
                (_sigma(Sigma::PHPH, th, ph + 1) -
                 _sigma(Sigma::PHPH, th, ph - 1)) /
                (2 * dPh);
        }
    }
}

void Conductance::_calcCoeff() {

    for (size_t th = 0; th < _nTh; th++) {
        for (size_t ph = 0; ph < _nPh; ph++) {
            double cos = std::cos((*_coords)(Ang::TH, th, ph));
            double sin = std::sin((*_coords)(Ang::TH, th, ph));
            double sin2 = std::pow(sin, 2);

            (*_coefficients)(Coeff::THTH, th, ph) =
                sin2 * _sigma(Sigma::THTH, th, ph);
            (*_coefficients)(Coeff::PHPH, th, ph) = _sigma(Sigma::PHPH, th, ph);
            (*_coefficients)(Coeff::TH, th, ph) =
                sin2 * _dSigma(DSigma::DTHTH_TH, th, ph) +
                sin * cos * _sigma(Sigma::THTH, th, ph) -
                sin * _dSigma(DSigma::DTHPH_PH, th, ph);
            (*_coefficients)(Coeff::PH, th, ph) =
                sin * _dSigma(DSigma::DTHPH_TH, th, ph) +
                _dSigma(DSigma::DPHPH_PH, th, ph);
        }
    }
}
