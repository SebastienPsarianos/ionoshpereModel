#include "Teuchos_RCPDecl.hpp"
#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <string>

class LegacyMHDConversion {
  public:
    static void processLegacyOutput(Teuchos::RCP<Coordinates>& coordinates,
                                    Teuchos::RCP<DipoleModel> dipoleModel,
                                    Teuchos::RCP<SolarModel> solarModel,
                                    Ionosphere::VectorRCP& sourceTerm,
                                    Ionosphere::MapRCP map,
                                    Ionosphere::CommRCP comm, int nTh, int nPh,
                                    const std::string& filename);
    static void getGridSize(std::string filename, int* nTh, int* nPh,
                            Ionosphere::CommRCP comm);
};
