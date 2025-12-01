#include "Grid.hpp"
#include <memory>
#include <string>

class LegacyMHDConversion {
  public:
    static int processLegacyOutput(std::string filename,
                                   std::shared_ptr<GridSet<Ang>> coordValues,
                                   std::shared_ptr<Grid> radCurrent, int nTh,
                                   int nPh);
    static int getGridSize(std::string filename, int* nTh, int* nPh);
};
