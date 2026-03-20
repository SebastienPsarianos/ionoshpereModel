#pragma once
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/conductance/EuvConductance.hpp"
#include "ionosphere/conductance/HardyConductance.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_RCP.hpp>
#include <nlohmann/json>

class EmpConductance {

  public:
    EmpConductance(Teuchos::RCP<Coordinates> coords, Ionosphere::MapRCP map,
                   Ionosphere::Scalar sig0);

    std::tuple<Ionosphere::MultiVectorRCP, Ionosphere::MultiVectorRCP,
               Ionosphere::MultiVectorRCP>
    computeConductance(int kp, double f107);

  private:
    Ionosphere::MultiVectorRCP
    _computeHppConductance(Ionosphere::MultiVectorRCP auroralConductance,
                           Ionosphere::MultiVectorRCP euvConductance);

    Ionosphere::MapRCP _map;
    Teuchos::RCP<Coordinates> _coords;

    HardyConductance hardyConductanceModel;
    EuvConductance euvConductanceModel;

    double _sig0;
};
