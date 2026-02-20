#pragma once
#include "ionosphere/utils/Grid.hpp"

#include <Eigen/Dense>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Tpetra_Map_fwd.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <nlohmann/json>

class EmpConductance {
  public:
    EmpConductance(MultiVectorRcp& sigma, MultiVectorRcp& coords, size_t nTh,
                   size_t nPh, double sig0, int year, int day, int month,
                   int hour, MapRcp map);

    void computeConductance(int kp, double f107);

  private:
    void _calcSigmaDer(Grid<DSigma>& dsigma, Grid<Sigma>& sigma,
                       Grid<GeoSph>& coords);

    void _computeAuroralConductance(int kp);
    void _computeEuvConductance(double f107);
    void _computeHppConductance();
    void _readAndSyncJson(Teuchos::RCP<const Teuchos::Comm<int>> comm);

    MapRcp _map;
    MultiVectorRcp _sigma;
    MultiVectorRcp _coords;
    MultiVectorRcp _euvConductance;
    MultiVectorRcp _auroralConductance;

    nlohmann::json _coefficientJson;

    size_t _nTh;
    size_t _nPh;

    double _utTime;
    double _ttTime;
    double _sig0;

    Eigen::Matrix3d _rotationMatrix;
};
