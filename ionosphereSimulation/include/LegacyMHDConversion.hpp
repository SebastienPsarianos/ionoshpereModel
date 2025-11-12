#include "Grid.hpp"
#include <string>

class LegacyMHDConversion {
  public:
    static int processLegacyOutput(std::string filename,
                                   GridSet<Ang>& coordValues, Grid& radCurrent,
                                   int nTh, int nPh);
    static int getGridSize(std::string filename, int* nTh, int* nPh);
};
