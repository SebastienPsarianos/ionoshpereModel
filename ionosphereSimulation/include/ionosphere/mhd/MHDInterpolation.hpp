#pragma once

#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_RCPDecl.hpp>

class MHDInterpolation {
  public:
    static void processLegacyOutput(Teuchos::RCP<Coordinates>& coordinates,
                                    Teuchos::RCP<DipoleModel> dipoleModel,
                                    Teuchos::RCP<SolarModel> solarModel,
                                    Ionosphere::VectorRCP& sourceTerm,
                                    Ionosphere::MapRCP map,
                                    Ionosphere::CommRCP comm, int nTh, int nPh,
                                    const Teuchos::ParameterList& ioParams);

    static void getGridSize(const Teuchos::ParameterList& ioParams, int* nTh,
                            int* nPh, Ionosphere::CommRCP comm);
};
