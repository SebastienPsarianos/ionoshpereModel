#include "ionosphere/utils/LegacyMHDConversion.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

int LegacyMHDConversion::processLegacyOutput(std::string filename,
                                             Grid<GeoSph>& coords,
                                             Grid<double>& radCurrent, int nTh,
                                             int nPh) {
    std::fstream jrData = std::fstream(filename);

    if (!jrData.is_open()) {
        std::cerr << "Error: Failed to open input file" << std::endl;
        return -1;
    }

    std::string line;
    std::getline(jrData, line); // Ignore the line with the grid sizes

    for (int th = 0; th < nTh; th++) {
        for (int ph = 0; ph < nPh; ph++) {
            if (std::getline(jrData, line)) {

                std::sscanf(line.c_str(), "%lf %lf %lf", &coords(th, ph).theta,
                            &coords(th, ph).phi, &radCurrent(th, ph));
            } else {
                throw std::length_error(
                    "File length doesn't match provided coordinates");
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
        std::sscanf(line.c_str(), "nTh: %d, nPh: %d", nTh, nPh);
    }
    return 0;
}
