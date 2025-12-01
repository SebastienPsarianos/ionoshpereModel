#pragma once
#include "Grid.hpp"
#include <memory>

class Conductance {
  public:
    Conductance(size_t nTh, size_t nPh, double sig0, double sigP, double sigH,
                std::shared_ptr<GridSet<Ang>> coords)
        : _coords(coords), _sigma(nTh, nPh),
          _hppSigma(nTh, nPh,
                    std::map<HppSigma, double>{{HppSigma::HALL, sigH},
                                               {HppSigma::PEDERSON, sigP},
                                               {HppSigma::PARALLEL, sig0}}),
          _dSigma(nTh, nPh) {
        _nTh = coords->nTh;
        _nPh = coords->nPh;
        _coefficients = std::make_unique<GridSet<Coeff>>(nTh, nPh);
    }
    std::shared_ptr<GridSet<Coeff>> calculateCoefficients();

  private:
    void _calcSigma();
    void _calcSigmaDer();
    void _calcCoeff();

    size_t _nTh;
    size_t _nPh;

    std::shared_ptr<GridSet<Ang>> _coords;
    GridSet<Sigma> _sigma;
    GridSet<HppSigma> _hppSigma;
    GridSet<DSigma> _dSigma;
    std::shared_ptr<GridSet<Coeff>> _coefficients;
};
