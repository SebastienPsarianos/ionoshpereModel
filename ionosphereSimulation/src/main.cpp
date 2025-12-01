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

    int nTh = 0;
    int nPh = 0;

    LegacyMHDConversion::getGridSize(argv[1], &nTh, &nPh);

    std::cout << "Sucessfully able to open radial current data, (" << nTh
              << ", " << nPh << ") grid detected" << std::endl;

    std::shared_ptr<GridSet<Ang>> coords =
        std::make_shared<GridSet<Ang>>(nTh, nPh);
    std::shared_ptr<Grid> radCurrent = std::make_shared<Grid>(nTh, nPh, 0.0);

    LegacyMHDConversion::processLegacyOutput(std::string(argv[1]), coords,
                                             radCurrent, nTh, nPh);

    std::cout << "PRINTING" << std::endl;
    std::ofstream out("../data/new_radCurrent.txt");
    if (out.is_open()) {
        std::cout << " OPEN" << std::endl;
        radCurrent->printWithCoords(out, *coords);
    }

    return 0;

    Conductance conductance = Conductance(nTh, nPh, SIG0, SIGP, SIGH, coords);
    std::shared_ptr<GridSet<Coeff>> kappa = conductance.calculateCoefficients();

    Solver solver = Solver(nTh, nPh, kappa, coords, radCurrent, GAUSS_SEIDEL);
    std::shared_ptr<Grid> potential = solver.calculatePotential();
    std::shared_ptr<GridSet<Ang>> eField = calculateEField(potential, coords);

    auto potentialOutput = std::ofstream("../data/solvedPotential.txt", 'w');
    auto eFieldOutput = std::ofstream("../data/solvedEField.txt", 'w');

    if (potentialOutput.is_open()) {
        potential->printWithCoords(potentialOutput, *coords);
    } else {
        std::cerr << "Unable to open potential file" << std::endl;
    }

    if (eFieldOutput.is_open()) {
        eField->printWithCoords(eFieldOutput, *coords);
    } else {
        std::cerr << "Unable to open eField file" << std::endl;
    }
}
