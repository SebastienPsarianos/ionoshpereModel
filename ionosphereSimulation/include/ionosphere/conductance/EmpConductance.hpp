#pragma once
#include "ionosphere/IonosphereTypes.hpp"

#include <Eigen/Dense>
#include <nlohmann/json>

class EmpConductance {

  public:
    EmpConductance(Ionosphere::MultiVectorRCP sigma,
                   Ionosphere::MultiVectorRCP coords, Ionosphere::MapRCP map,
                   Ionosphere::Scalar sig0, size_t nTh, size_t nPh, int year,
                   int day, int month, int hour);

    void computeConductance(int kp, double f107);

  private:
    void _computeAuroralConductance(int kp);
    void _computeEuvConductance(double f107);
    void _computeHppConductance();
    void _readAndSyncJson(Teuchos::RCP<const Teuchos::Comm<int>> comm);

    Ionosphere::MapRCP _map;
    Ionosphere::MultiVectorRCP _sigma;
    Ionosphere::MultiVectorRCP _coords;
    Ionosphere::MultiVectorRCP _euvConductance;
    Ionosphere::MultiVectorRCP _auroralConductance;

    nlohmann::json _coefficientJson;

    size_t _nTh;
    size_t _nPh;

    double _utTime;
    double _ttTime;
    double _sig0;

    Eigen::Matrix3d _rotationMatrix;
};
