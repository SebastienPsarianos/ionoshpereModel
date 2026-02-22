#include "ionosphere/solver/TlSolver.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include <BelosLinearProblem.hpp> // Good practice to include this explicitly too
#include <BelosTpetraAdapter.hpp>

#include "Teuchos_DataAccess.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_Operator.hpp"
#include "ionosphere/utils/Grid.hpp"
#include <cmath>

// TODO: Figure out if this needs to be a class
TlSolver::TlSolver(size_t nTh, size_t nPh, double dPh, double dTh, MapRcp map)
    : Solver(nTh, nPh), _dTh(dTh), _dPh(dPh), _map(map) {}

VectorRcp TlSolver::calculatePotential(MultiVectorRcp conductance,
                                       MultiVectorRcp coords,
                                       VectorRcp sourceTerm) {
    MultiVectorRcp coefficients = _calculateCoefficients(conductance, coords);
    CrsMatrixRcp A = _buildGrid(coords, coefficients);
    auto x = Teuchos::rcp(
        new Tpetra::Vector<double, int, long long>(*sourceTerm, Teuchos::Copy));

    x->putScalar(0);
    auto problem = Teuchos::rcp(
        new Belos::LinearProblem<double,
                                 Tpetra::MultiVector<double, int, long long>,
                                 Tpetra::Operator<double, int, long long>>(
            A, x, sourceTerm));
    if (!problem->setProblem()) {
        std::cout << "EEEEE";
    }

    auto solverParams = Teuchos::parameterList();
    solverParams->set("Maximum Iterations", 1000);
    solverParams->set("Convergence Tolerance", 1e-8);

    Belos::BlockGmresSolMgr<double, Tpetra::MultiVector<double, int, long long>,
                            Tpetra::Operator<double, int, long long>>
        solver(problem, solverParams);

    Belos::ReturnType result = solver.solve();

    if (result == Belos::Converged) {
        std::cout << "WOHOOOO";
    } else {
        std::cout << "AAAAAAAAA";
    }

    return x;
}

CrsMatrixRcp TlSolver::_buildGrid(MultiVectorRcp coords,
                                  MultiVectorRcp coefficients) {

    auto ththCoefficients = coefficients->getDataNonConst(0);
    auto phphCoefficients = coefficients->getDataNonConst(1);
    auto thCoefficients = coefficients->getDataNonConst(2);
    auto phCoefficients = coefficients->getDataNonConst(3);

    auto th = coords->getDataNonConst(0);
    auto ph = coords->getDataNonConst(1);

    auto A = Teuchos::rcp(new Tpetra::CrsMatrix<>(_map, 5));
    auto myGridPoints = _map->getMyGlobalIndices();

    for (size_t i = 0; i < myGridPoints.size(); i++) {
        auto gridPoint = myGridPoints[i];

        size_t theta = gridPoint % _nTh;
        size_t phi = gridPoint / _nTh;

        Teuchos::Array<long long> matrixIdcs;
        Teuchos::Array<double> vals;

        matrixIdcs.push_back(gridPoint);
        vals.push_back(-2 * ththCoefficients[i] / (_dTh * _dTh) -
                       2 * phphCoefficients[i] / (_dPh * _dPh));

        size_t left = gridPoint - _nTh;
        size_t right = gridPoint + _nTh;
        size_t up = gridPoint - 1;
        size_t down = gridPoint + 1;

        if (theta == 0) {
            size_t oppositePhi = (phi + _nPh / 2) % _nPh;
            up = oppositePhi * _nTh;
        } else if (theta == _nTh - 1) {
            size_t oppositePhi = (phi + _nPh / 2) % _nPh;
            down = oppositePhi * _nTh + theta;
        }

        if (phi == 0) {
            left = theta + (_nPh - 1) * _nTh;
        } else if (phi == _nPh - 1) {
            right = theta;
        }

        // Up
        matrixIdcs.push_back(up);
        vals.push_back(ththCoefficients[i] / (_dTh * _dTh) -
                       thCoefficients[i] / (2 * _dTh));

        // Down
        matrixIdcs.push_back(down);
        vals.push_back(ththCoefficients[i] / (_dTh * _dTh) +
                       thCoefficients[i] / (2 * _dTh));

        // Left
        matrixIdcs.push_back(left);
        vals.push_back(phphCoefficients[i] / (_dPh * _dPh) -
                       phCoefficients[i] / (2 * _dPh));

        // Right
        matrixIdcs.push_back(right);
        vals.push_back(phphCoefficients[i] / (_dPh * _dPh) +
                       phCoefficients[i] / (2 * _dPh));

        A->insertGlobalValues(gridPoint, matrixIdcs, vals);
    }
    A->fillComplete();

    return A;
}

MultiVectorRcp TlSolver::_calculateCoefficients(MultiVectorRcp conductance,
                                                MultiVectorRcp coords) {
    auto D_th = Teuchos::rcp(new Tpetra::CrsMatrix<>(_map, 3));
    auto D_ph = Teuchos::rcp(new Tpetra::CrsMatrix<>(_map, 3));
    auto myGridPointsGlobal = _map->getMyGlobalIndices();

    // Build first derivative matrices
    for (size_t i = 0; i < myGridPointsGlobal.size(); i++) {
        long long currentPoint = myGridPointsGlobal[i];

        size_t theta = currentPoint % _nTh;
        size_t phi = currentPoint / _nTh;

        Teuchos::Array<long long> matrixIdcs_th;
        Teuchos::Array<double> vals_th;
        Teuchos::Array<long long> matrixIdcs_ph;
        Teuchos::Array<double> vals_ph;

        size_t left = currentPoint - _nTh;
        size_t right = currentPoint + _nTh;
        size_t up = currentPoint - 1;
        size_t down = currentPoint + 1;

        if (theta == 0) {
            size_t oppositePhi = (phi + _nPh / 2) % _nPh;
            up = oppositePhi * _nTh;
        } else if (theta == _nTh - 1) {
            size_t oppositePhi = (phi + _nPh / 2) % _nPh;
            down = oppositePhi * _nTh + theta;
        }

        if (phi == 0) {
            left = theta + (_nPh - 1) * _nTh;
        } else if (phi == _nPh - 1) {
            right = theta;
        }

        matrixIdcs_th.push_back(down);
        vals_th.push_back(1.0 / (2 * _dTh));
        matrixIdcs_th.push_back(up);
        vals_th.push_back(-1.0 / (2 * _dTh));

        matrixIdcs_ph.push_back(right);
        vals_ph.push_back(1.0 / (2 * _dPh));
        matrixIdcs_ph.push_back(left);
        vals_ph.push_back(-1.0 / (2 * _dPh));

        D_th->insertGlobalValues(currentPoint, matrixIdcs_th, vals_th);
        D_ph->insertGlobalValues(currentPoint, matrixIdcs_ph, vals_ph);
    }
    D_th->fillComplete();
    D_ph->fillComplete();

    // Create multivectors containing the terms in the coefficients where the
    // derivative is taken. Stored in the order given in O.Amm equation 2
    auto thDerToTake =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(_map, 2));
    auto phDerToTake =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(_map, 2));

    auto thd1 = thDerToTake->getDataNonConst(0);
    auto thd2 = thDerToTake->getDataNonConst(1);
    auto phd1 = phDerToTake->getDataNonConst(0);
    auto phd2 = phDerToTake->getDataNonConst(1);

    auto thVals = coords->getDataNonConst(0);
    auto phVals = coords->getDataNonConst(1);
    auto pedersonVals = conductance->getDataNonConst(0);
    auto hallVals = conductance->getDataNonConst(1);
    auto parallelVals = conductance->getDataNonConst(2);

    for (int i = 0; i < thVals.size(); i++) {
        // Eqn 21, 22 from goodman
        double cos_e =
            (-2 * std::cos(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));
        double sin_e =
            (std::sin(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));

        double C =
            parallelVals[i] * cos_e * cos_e + pedersonVals[i] * sin_e * sin_e;

        thd1[i] = parallelVals[i] * pedersonVals[i] / C;
        thd2[i] = (-parallelVals[i] * hallVals[i] * cos_e) /
                  (C * std::sin(thVals[i]));
        phd1[i] = parallelVals[i] * hallVals[i] * cos_e / C;
        phd2[i] =
            pedersonVals[i] + (hallVals[i] * hallVals[i] * sin_e * sin_e) / C;
    }

    // Apply derivatives to terms as mentioned above
    auto thDer =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(_map, 2));
    auto phDer =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(_map, 2));

    D_th->apply(*thDerToTake, *thDer);
    D_ph->apply(*phDerToTake, *phDer);

    // Build PDE-coefficient multi-vector

    auto coefficients =
        Teuchos::rcp(new Tpetra::MultiVector<double, int, long long>(_map, 4));

    auto kappa_thth = coefficients->getDataNonConst(0);
    auto kappa_phph = coefficients->getDataNonConst(1);
    auto kappa_th = coefficients->getDataNonConst(2);
    auto kappa_ph = coefficients->getDataNonConst(3);

    thd1 = thDer->getDataNonConst(0);
    thd2 = thDer->getDataNonConst(1);
    phd1 = phDer->getDataNonConst(0);
    phd2 = phDer->getDataNonConst(1);

    for (int i = 0; i < kappa_thth.size(); i++) {
        // Eqn 21, 22 from goodman
        double cos_e =
            (-2 * std::cos(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));
        double sin_e =
            (std::sin(thVals[i])) /
            std::sqrt((1 + 3 * std::cos(thVals[i]) * std::cos(thVals[i])));

        double C =
            parallelVals[i] * cos_e * cos_e + pedersonVals[i] * sin_e * sin_e;

        double sin = std::sin(thVals[i]);
        double cot = 1.0 / std::tan(thVals[i]);

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
