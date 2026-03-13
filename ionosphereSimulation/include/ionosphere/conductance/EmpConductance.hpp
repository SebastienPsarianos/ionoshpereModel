#pragma once
#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_RCP.hpp>
#include <nlohmann/json>

class EmpConductance {

  public:
    EmpConductance(Teuchos::RCP<Coordinates> coords, Ionosphere::MapRCP map,
                   Ionosphere::Scalar sig0);

    Ionosphere::MultiVectorRCP computeConductance(int kp, double f107);

  private:
    void _computeAuroralConductance(int kp);
    void _computeEuvConductance(double f107);
    void _computeHppConductance();
    void _readAndSyncJson(Ionosphere::CommRCP comm);

    Ionosphere::MapRCP _map;
    Ionosphere::MultiVectorRCP _sigma;
    Ionosphere::MultiVectorRCP _euvConductance;
    Ionosphere::MultiVectorRCP _auroralConductance;
    Ionosphere::MultiVectorRCP _conductance;

    nlohmann::json _coefficientJson;

    double _sig0;

    Teuchos::RCP<Coordinates> _coords;
};
