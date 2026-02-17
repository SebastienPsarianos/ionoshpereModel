#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/postProcessing/EField.hpp"
#include "ionosphere/solver/GsSolver.hpp"
#include "ionosphere/utils/Coordinates.hpp"
#include "ionosphere/utils/Grid.hpp"
#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <fstream>
#include <iostream>

double SIG0 = 20.00;
double SIGP = 2.00;
double SIGH = 2.00;

int YEAR = 2026;
int MONTH = 2;
int DAY = 5;
int HOUR = 5;

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
    Grid<GeoSph> coords = Grid<GeoSph>(nTh, nPh);
    Grid<double> radCurrent = Grid<double>(nTh, nPh, 0.0);
    Grid<Sigma> conductance = Grid<Sigma>(nTh, nPh);
    Grid<DSigma> dConductance = Grid<DSigma>(nTh, nPh);
    Grid<double> potential = Grid<double>(nTh, nPh, 0.0);
    Grid<GeoSph> eField = Grid<GeoSph>(nTh, nPh);

    // TODO: Coupling procedure
    int processRes = LegacyMHDConversion::processLegacyOutput(
        std::string(argv[1]), coords, radCurrent, nTh, nPh);

    if (processRes > 0) {
        std::cerr << "Failed to process simulation input data" << std::endl;
        return 1;
    }

    EmpConductance(nTh, nPh, SIG0, 2026, 2, 5, 5)
        .computeConductance(conductance, dConductance, coords);

    // TODO: Implement trilinos solver
    GsSolver(potential, nTh, nPh, coords, radCurrent, conductance, dConductance)
        .calculatePotential();

    // TODO: Proper post-processing, including J and others
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
