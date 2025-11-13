#include "Conductance.hpp"
#include "Coordinates.hpp"
#include "Grid.hpp"
#include "LegacyMHDConversion.hpp"
#include "Solver.hpp"
#include <iostream>

double SIG0 = 20.00;
double SIGP = 2.00;
double SIGH = 2.00;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Missing filename for input data" << std::endl;
        return -1;
    }

    int nTh = 0;
    int nPh = 0;

    LegacyMHDConversion::getGridSize(argv[1], &nTh, &nPh);

    std::shared_ptr<GridSet<Ang>> coords =
        std::make_shared<GridSet<Ang>>(nTh, nPh);
    std::shared_ptr<Grid> radCurrent = std::make_shared<Grid>(nTh, nPh, 0.0);

    LegacyMHDConversion::processLegacyOutput(std::string(argv[1]), coords,
                                             radCurrent, nTh, nPh);

    Conductance conductance = Conductance(nTh, nPh, SIG0, SIGP, SIGH, coords);
    std::shared_ptr<GridSet<Coeff>> kappa = conductance.calculateCoefficients();

    Solver solver = Solver(nTh, nPh, kappa, coords, radCurrent, GAUSS_SEIDEL);
    std::shared_ptr<Grid> potential = solver.calculatePotential();

    std::cout << potential;
}
