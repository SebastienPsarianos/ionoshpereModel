#pragma once
#include "Teuchos_ParameterList.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/conductance/EuvConductance.hpp"
#include "ionosphere/conductance/HardyConductance.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_RCP.hpp>
#include <nlohmann/json>

class EmpConductance {

  public:
    EmpConductance(Teuchos::RCP<Coordinates> coords, Ionosphere::MapRCP map,
                   const Teuchos::ParameterList& conductanceParams);

    std::tuple<Ionosphere::MultiVectorRCP, Ionosphere::MultiVectorRCP,
               Ionosphere::MultiVectorRCP>
    computeConductance();

  private:
    Ionosphere::MultiVectorRCP
    _computeHppConductance(Ionosphere::MultiVectorRCP auroralConductance,
                           Ionosphere::MultiVectorRCP euvConductance);

    Ionosphere::MapRCP _map;
    Teuchos::RCP<Coordinates> _coords;

    HardyConductance hardyConductanceModel;
    EuvConductance euvConductanceModel;

    int _kp;
    double _f107;
    double _sig0;
};
