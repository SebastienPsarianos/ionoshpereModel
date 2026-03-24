#include "ionosphere/solver/OAmmEqn.hpp"
#include "ionosphere/TrilinosAliases.hpp"

using namespace Ionosphere;

OAmmEqn::OAmmEqn(const Teuchos::ParameterList& equationParams,
                 Teuchos::RCP<Coordinates> coords, MultiVectorRCP conductance,
                 VectorRCP radCurrent, MapRCP map)
    : _coords(coords), _conductance(conductance), _radCurrent(radCurrent),
      _map(map),
      _radiusIonosphere2((equationParams.get<double>("ionosphere_altitude_m") +
                          equationParams.get<double>("earth_radius_m")) *
                         (equationParams.get<double>("ionosphere_altitude_m") +
                          equationParams.get<double>("earth_radius_m"))) {}

// TODO: Got to verify all of this
MatrixRCP OAmmEqn::assembleMatrix() {
    using std::sin;

    Ionosphere::MultiVectorRCP coefficients = _calculateCoefficients();

    auto ththCoefficients = coefficients->getDataNonConst(0);
    auto phphCoefficients = coefficients->getDataNonConst(1);
    auto thCoefficients = coefficients->getDataNonConst(2);
    auto phCoefficients = coefficients->getDataNonConst(3);

    auto pedersonVals = _conductance->getDataNonConst(0);
    auto hallVals = _conductance->getDataNonConst(1);
    auto parallelVals = _conductance->getDataNonConst(2);

    auto coordsMv = _coords->multiVector();
    auto th = coordsMv->getDataNonConst(0);
    auto ph = coordsMv->getDataNonConst(1);

    auto A = Teuchos::rcp(
        new Matrix(_map, std::max(size_t(5), size_t(_coords->nPh + 1))));
    auto myGridPoints = _map->getMyGlobalIndices();

    // Set phi == 0 on the equator to be a gauge condition
    for (GlobalOrd i = 0; i < static_cast<GlobalOrd>(myGridPoints.size());
         i++) {
        GlobalOrd gridPoint = myGridPoints[i];
        auto [theta, phi] = _coords->globalIdx2ThetaPhi(gridPoint);

        Teuchos::Array<GlobalOrd> matrixIdcs;
        Teuchos::Array<Scalar> vals;

        // Pin the gauge point on the equator and then skip
        if (theta == _coords->nTh / 2 && phi == 0) {
            matrixIdcs.push_back(gridPoint);
            vals.push_back(1);
            A->insertGlobalValues(gridPoint, matrixIdcs, vals);
            continue;
        }

        // Applying the polar cap flux boundary condition (see notes)
        if (theta != 0 && theta != _coords->nTh - 1) {
            matrixIdcs.push_back(gridPoint);
            vals.push_back(
                -2 * ththCoefficients[i] / (_coords->dTh * _coords->dTh) -
                2 * phphCoefficients[i] / (_coords->dPh * _coords->dPh));

            GlobalOrd left = gridPoint - _coords->nTh;
            GlobalOrd right = gridPoint + _coords->nTh;

            GlobalOrd up = gridPoint - 1;
            GlobalOrd down = gridPoint + 1;

            if (phi == 0) {
                left = theta + (_coords->nPh - 1) * _coords->nTh;
            } else if (phi == _coords->nPh - 1) {
                right = theta;
            }

            // Up
            matrixIdcs.push_back(up);
            vals.push_back(ththCoefficients[i] / (_coords->dTh * _coords->dTh) -
                           thCoefficients[i] / (2 * _coords->dTh));

            // Down
            matrixIdcs.push_back(down);
            vals.push_back(ththCoefficients[i] / (_coords->dTh * _coords->dTh) +
                           thCoefficients[i] / (2 * _coords->dTh));

            // Left
            matrixIdcs.push_back(left);
            vals.push_back(phphCoefficients[i] / (_coords->dPh * _coords->dPh) -
                           phCoefficients[i] / (2 * _coords->dPh));

            // Right
            matrixIdcs.push_back(right);
            vals.push_back(phphCoefficients[i] / (_coords->dPh * _coords->dPh) +
                           phCoefficients[i] / (2 * _coords->dPh));

            A->insertGlobalValues(gridPoint, matrixIdcs, vals);
        }
    }

    return A;
}

VectorRCP OAmmEqn::assembleRHS() {

    auto rhs = Teuchos::rcp(new Vector(*_radCurrent, Teuchos::Copy));
    rhs->scale(_radiusIonosphere2);

    return rhs;
}

VectorRCP OAmmEqn::initialGuess() {
    auto x = Teuchos::rcp(new Vector(_map));
    x->putScalar(0.0);
    return x;
}

MultiVectorRCP OAmmEqn::_calculateCoefficients() {
    auto D_th = Teuchos::rcp(new Matrix(_map, 4));
    auto D_ph = Teuchos::rcp(new Matrix(_map, 3));

    auto myGridPointsGlobal = _map->getMyGlobalIndices();
    // Build first derivative matrices
    for (LocalOrd i = 0; i < static_cast<LocalOrd>(myGridPointsGlobal.size());
         i++) {
        GlobalOrd currentGid = myGridPointsGlobal[i];
        auto [theta, phi] = _coords->globalIdx2ThetaPhi(currentGid);

        Teuchos::Array<GlobalOrd> matrixIdcs_th;
        Teuchos::Array<Scalar> vals_th;
        Teuchos::Array<GlobalOrd> matrixIdcs_ph;
        Teuchos::Array<Scalar> vals_ph;

        GlobalOrd left = currentGid - _coords->nTh;
        GlobalOrd right = currentGid + _coords->nTh;
        GlobalOrd up = currentGid - 1;
        GlobalOrd down = currentGid + 1;

        if (phi == 0) {
            left = theta + (_coords->nPh - 1) * _coords->nTh;
        } else if (phi == _coords->nPh - 1) {
            right = theta;
        }

        if (theta == 0 || theta == _coords->nTh - 1) {
            matrixIdcs_th.push_back(currentGid);
            vals_th.push_back(0.0);
        } else if (theta == 1) {
            // One sided stencil for north pole (second order)
            GlobalOrd pt1 = phi * _coords->nTh + 1;
            GlobalOrd pt2 = phi * _coords->nTh + 2;
            GlobalOrd pt3 = phi * _coords->nTh + 3;
            matrixIdcs_th.push_back(pt1);
            vals_th.push_back(-3.0 / (2.0 * _coords->dTh));
            matrixIdcs_th.push_back(pt2);
            vals_th.push_back(4.0 / (2.0 * _coords->dTh));
            matrixIdcs_th.push_back(pt3);
            vals_th.push_back(-1.0 / (2.0 * _coords->dTh));
        } else if (theta == _coords->nTh - 2) {
            // One sided stencil for south pole (second order)
            GlobalOrd pt1 = phi * _coords->nTh + theta;
            GlobalOrd pt2 = phi * _coords->nTh + theta - 1;
            GlobalOrd pt3 = phi * _coords->nTh + theta - 2;
            matrixIdcs_th.push_back(pt1);
            vals_th.push_back(3.0 / (2.0 * _coords->dTh));
            matrixIdcs_th.push_back(pt2);
            vals_th.push_back(-4.0 / (2.0 * _coords->dTh));
            matrixIdcs_th.push_back(pt3);
            vals_th.push_back(1.0 / (2.0 * _coords->dTh));
        } else {
            // Standard Derrivative stencil for everywhere else
            matrixIdcs_th.push_back(down);
            vals_th.push_back(1.0 / (2.0 * _coords->dTh));
            matrixIdcs_th.push_back(up);
            vals_th.push_back(-1.0 / (2.0 * _coords->dTh));
        }

        if (theta == 0 || theta == _coords->nTh - 1) {
            matrixIdcs_ph.push_back(currentGid);
            vals_ph.push_back(0.0);
        } else {
            matrixIdcs_ph.push_back(right);
            vals_ph.push_back(1.0 / (2 * _coords->dPh));
            matrixIdcs_ph.push_back(left);
            vals_ph.push_back(-1.0 / (2 * _coords->dPh));
        }
        D_th->insertGlobalValues(currentGid, matrixIdcs_th, vals_th);
        D_ph->insertGlobalValues(currentGid, matrixIdcs_ph, vals_ph);
    }
    D_th->fillComplete();
    D_ph->fillComplete();

    // Create multivectors containing the terms in the coefficients where the
    // derivative is taken. Stored in the order given in O.Amm equation 2
    auto thDerToTake = Teuchos::rcp(new MultiVector(_map, 2));
    auto phDerToTake = Teuchos::rcp(new MultiVector(_map, 2));

    auto thd1 = thDerToTake->getDataNonConst(0);
    auto thd2 = thDerToTake->getDataNonConst(1);
    auto phd1 = phDerToTake->getDataNonConst(0);
    auto phd2 = phDerToTake->getDataNonConst(1);

    auto coordsMv = _coords->multiVector();
    auto thVals = coordsMv->getDataNonConst(0);
    auto phVals = coordsMv->getDataNonConst(1);

    auto pedersonVals = _conductance->getDataNonConst(0);
    auto hallVals = _conductance->getDataNonConst(1);
    auto parallelVals = _conductance->getDataNonConst(2);

    for (LocalOrd i = 0; i < static_cast<LocalOrd>(thVals.size()); i++) {
        GlobalOrd thetaGlobal = myGridPointsGlobal[i] % _coords->nTh;

        if (thetaGlobal == 0 || thetaGlobal == _coords->nTh - 1) {
            thd1[i] = 0.0;
            thd2[i] = 0.0;
            phd1[i] = 0.0;
            phd2[i] = 0.0;
            continue;
        }

        // Eqn 21, 22 from goodman
        Scalar cos_e =
            (-2 * std::cos(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));
        Scalar sin_e =
            (std::sin(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));

        Scalar C =
            parallelVals[i] * cos_e * cos_e + pedersonVals[i] * sin_e * sin_e;

        thd1[i] = parallelVals[i] * pedersonVals[i] / C;
        thd2[i] = (-parallelVals[i] * hallVals[i] * cos_e) /
                  (C * std::sin(thVals[i]));
        phd1[i] = parallelVals[i] * hallVals[i] * cos_e / C;
        phd2[i] =
            pedersonVals[i] + (hallVals[i] * hallVals[i] * sin_e * sin_e) / C;
    }

    // Apply derivatives to terms as mentioned above
    auto thDer = Teuchos::rcp(new MultiVector(_map, 2));
    auto phDer = Teuchos::rcp(new MultiVector(_map, 2));

    D_th->apply(*thDerToTake, *thDer);
    D_ph->apply(*phDerToTake, *phDer);

    // Build PDE-coefficient multi-vector
    auto coefficients = Teuchos::rcp(new MultiVector(_map, 4));

    auto kappa_thth = coefficients->getDataNonConst(0);
    auto kappa_phph = coefficients->getDataNonConst(1);
    auto kappa_th = coefficients->getDataNonConst(2);
    auto kappa_ph = coefficients->getDataNonConst(3);

    thd1 = thDer->getDataNonConst(0);
    thd2 = thDer->getDataNonConst(1);
    phd1 = phDer->getDataNonConst(0);
    phd2 = phDer->getDataNonConst(1);

    for (LocalOrd i = 0; i < static_cast<LocalOrd>(kappa_thth.size()); i++) {
        // Eqn 21, 22 from goodman
        GlobalOrd thetaGlobal = myGridPointsGlobal[i] % _coords->nTh;

        if (thetaGlobal == 0 || thetaGlobal == _coords->nTh - 1) {
            kappa_thth[i] = 0.0;
            kappa_phph[i] = 0.0;
            kappa_th[i] = 0.0;
            kappa_ph[i] = 0.0;
            continue;
        }
        Scalar cos_e =
            (-2 * std::cos(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));
        Scalar sin_e =
            (std::sin(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));

        Scalar C =
            parallelVals[i] * cos_e * cos_e + pedersonVals[i] * sin_e * sin_e;

        Scalar sin = std::sin(thVals[i]);
        Scalar cot = 1.0 / std::tan(thVals[i]);

        kappa_thth[i] = parallelVals[i] * pedersonVals[i] / C;
        kappa_phph[i] = (pedersonVals[i] +
                         (hallVals[i] * hallVals[i] * sin_e * sin_e) / C) /
                        (sin * sin);
        kappa_th[i] = thd1[i] + (cot * parallelVals[i] * pedersonVals[i]) / C +
                      (1 / sin) * phd1[i];
        kappa_ph[i] = thd2[i] + (1 / (sin * sin)) * phd2[i] -
                      (parallelVals[i] * hallVals[i] * cos_e * cot) / (C * sin);
    }

    return coefficients;
}
