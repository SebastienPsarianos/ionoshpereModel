#include "Conductance.hpp"
#include "Coordinates.hpp"
#include "EField.hpp"
#include "Grid.hpp"
#include "LegacyMHDConversion.hpp"
#include "Solver.hpp"
#include <fstream>
#include <iostream>

double SIG0 = 20.00;
double SIGP = 2.00;
double SIGH = 2.00;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Missing filename for Radial Current data" << std::endl;
        return -1;
    }

    int nTh = -1;
    int nPh = -1;

    int gridSizeRes = LegacyMHDConversion::getGridSize(argv[1], &nTh, &nPh);
    if (gridSizeRes > 0 || nTh <= 0 || nPh <= 0) {
        std::cerr << "Error reading file ending execution" << std::endl;
        return 1;
    }
    std::cout << "Sucessfully able to open radial current data, (" << nTh
              << ", " << nPh << ") grid detected" << std::endl;

    // creating required grids
    Grid<ThPh> coords = Grid<ThPh>(nTh, nPh);
    Grid<double> radCurrent = Grid<double>(nTh, nPh, 0.0);
    Grid<Coeff> kappa = Grid<Coeff>(nTh, nPh);
    Grid<Sigma> conductance = Grid<Sigma>(nTh, nPh);
    Grid<double> potential = Grid<double>(nTh, nPh, 0.0);
    Grid<ThPh> eField = Grid<ThPh>(nTh, nPh);

    int processRes = LegacyMHDConversion::processLegacyOutput(
        std::string(argv[1]), coords, radCurrent, nTh, nPh);

    if (processRes > 0) {
        std::cerr << "Failed to process simulation input data" << std::endl;
        return 1;
    }

    // Calculate the coefficients
    Conductance(kappa, conductance, nTh, nPh, SIG0, SIGP, SIGH, coords)
        .calculateCoefficients();

    // TODO: Calculate the whole source term here so we don't need to
    // recalculate in gauss seidel

    // Calculate the potential
    Solver(potential, nTh, nPh, kappa, coords, radCurrent, conductance,
           GAUSS_SEIDEL)
        .calculatePotential();

    // Calculate the electric field
    calculateEField(eField, potential, coords);

    std::ofstream potentialOutput("../data/solvedPotential.txt");
    std::ofstream eFieldOutput("../data/solvedEField.txt");

    if (potentialOutput.is_open()) {
        potential.printWithCoords(potentialOutput, coords);
    } else {
        std::cerr << "Unable to open potential file" << std::endl;
    }

    if (eFieldOutput.is_open()) {
        eField.printWithCoords(eFieldOutput, coords);
    } else {
        std::cerr << "Unable to open eField file" << std::endl;
    }
}
