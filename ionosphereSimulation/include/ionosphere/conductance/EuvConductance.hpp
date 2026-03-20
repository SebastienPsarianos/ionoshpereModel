#pragma once
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_RCP.hpp>

class EuvConductance {

  public:
    EuvConductance(Teuchos::RCP<Coordinates> coords, Ionosphere::MapRCP map);

    Ionosphere::MultiVectorRCP computeEuvConductance(double f107);

  private:
    Ionosphere::MapRCP _map;
    Ionosphere::MultiVectorRCP _euvConductance;
    Teuchos::RCP<Coordinates> _coords;
};
