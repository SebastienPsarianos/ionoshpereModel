#include "Coordinates.hpp"
#include "Grid.hpp"
#include "LegacyMHDConversion.hpp"
#include <iostream>

std::map<HppSigma, double> HPPVALS = {{HppSigma::PARALLEL, 20.00},
                                      {HppSigma::HALL, 2.00},
                                      {HppSigma::PEDERSON, 2.00}};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Missing filename for input data" << std::endl;
        return -1;
    }

    int nTh = 0;
    int nPh = 0;

    LegacyMHDConversion::getGridSize(argv[1], &nTh, &nPh);

    GridSet<Sigma> sigma = GridSet<Sigma>(nTh, nPh);
    GridSet<HppSigma> hppSigma = GridSet<HppSigma>(nTh, nPh, HPPVALS);
    GridSet<DSigma> dSigma = GridSet<DSigma>(nTh, nPh);
    GridSet<Ang> coords = GridSet<Ang>(nTh, nPh);
    Grid radCurrent = Grid(nTh, nPh, 0.0);

    LegacyMHDConversion::processLegacyOutput(std::string(argv[1]), coords,
                                             radCurrent, nTh, nPh);
}
