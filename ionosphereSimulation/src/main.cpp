#include "Coordinates.hpp"
#include "Grid.hpp"
#include <iostream>

const size_t NTH = 10;
const size_t NPH = 20;
std::map<HppSigma, double> HPPVALS = {{HppSigma::PARALLEL, 20.00},
                                      {HppSigma::HALL, 2.00},
                                      {HppSigma::PEDERSON, 2.00}};

int main() {
    GridSet<Sigma> sigma = GridSet<Sigma>(NTH, NPH);
    GridSet<HppSigma> hppSigma = GridSet<HppSigma>(NTH, NPH, HPPVALS);
    GridSet<DSigma> dSigma = GridSet<DSigma>(NTH, NPH);
    GridSet<Coords> coords = GridSet<Coords>(NTH, NPH);

    Coordinates::calculateCoords(coords);

    std::cout << "THETA" << std::endl;

    for (size_t i = 0; i < NTH; i++) {
        for (size_t j = 0; j < NPH; j++) {
            std::cout << coords(Coords::TH, i, j) << " ";
        }
        std::cout << "\n";
    }

    std::cout << "PHI" << std::endl;
    for (size_t i = 0; i < NTH; i++) {
        for (size_t j = 0; j < NPH; j++) {
            std::cout << coords(Coords::PH, i, j) << " ";
        }
        std::cout << "\n";
    }
}
