#pragma once
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"

#include <Teuchos_RCP.hpp>
#include <nlohmann/json>

class HardyConductance {
  public:
    HardyConductance(Teuchos::RCP<Coordinates> coords, Ionosphere::MapRCP map);
    Ionosphere::MultiVectorRCP computeAuroralConductance(int kp);

  private:
    void _readAndSyncJson(Ionosphere::CommRCP comm);
    double _fourierSeries(nlohmann::json coefficients, double mlt);
    double _epsteinFunction(double h, double r, double h0, double S1,
                            double S2);

    Ionosphere::MapRCP _map;
    Ionosphere::MultiVectorRCP _auroralConductance;
    Teuchos::RCP<Coordinates> _coords;

    nlohmann::json _coefficientJson;
};
