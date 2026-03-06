#include "Teuchos_RCPDecl.hpp"
#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/utils/Coordinates.hpp"

#include <string>

class LegacyMHDConversion {
  public:
    static void processLegacyOutput(Teuchos::RCP<Coordinates>& coordinates,
                                    Ionosphere::VectorRCP& sourceTerm,
                                    Ionosphere::MapRCP map,
                                    Ionosphere::CommRCP comm, int nTh, int nPh,
                                    const std::string& filename);
    static void getGridSize(std::string filename, int* nTh, int* nPh,
                            Ionosphere::CommRCP comm);
};
