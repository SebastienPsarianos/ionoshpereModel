#include "LegacyMHDConversion.hpp"
#include <fstream>
#include <iostream>

int LegacyMHDConversion::processLegacyOutput(
    std::string filename, std::shared_ptr<GridSet<Ang>> coordValues,
    std::shared_ptr<Grid> radCurrent, int nTh, int nPh) {
    std::fstream dataFile = std::fstream(filename);

    if (!dataFile.is_open()) {
        std::cerr << "Error: Failed to open input file" << std::endl;
        return -1;
    }

    dataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string line;
    for (int th = 0; th < nTh; th++) {
        for (int ph = 0; ph < nPh; ph++) {
            if (!dataFile.eof()) {
                std::getline(dataFile, line);
                std::sscanf(line.c_str(), "%lf %lf %lf",
                            &(*coordValues)(Ang::TH, th, ph),
                            &(*coordValues)(Ang::PH, th, ph),
                            &(*radCurrent)(th, ph));

            } else {
                std::cerr << "Invalid Data File\n";
                return -1;
            }
        }
    }

    return 0;
}

int LegacyMHDConversion::getGridSize(std::string filename, int* nTh, int* nPh) {
    std::fstream dataFile = std::fstream(filename);
    if (!dataFile.is_open()) {
        std::cerr << "Error: Failed to open input file" << std::endl;
        return -1;
    }

    std::string line;
    if (std::getline(dataFile, line)) {
        std::sscanf(line.c_str(), "th: %d, ph: %d", nTh, nPh);
    }
    return 0;
}
