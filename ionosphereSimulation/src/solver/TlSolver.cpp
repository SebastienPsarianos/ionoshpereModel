#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/Constants.hpp"
#include "ionosphere/IonosphereTypes.hpp"

#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Tpetra_Operator.hpp>
#include <cmath>
#include <stdexcept>

using namespace Ionosphere;

// TODO: Figure out if this needs to be a class
TlSolver::TlSolver(Teuchos::RCP<Coordinates> coords, MapRCP map)
    : Solver(coords->nTh, coords->nPh), _dTh(coords->dTh), _dPh(coords->dPh),
      _coords(coords), _map(map) {}

VectorRCP TlSolver::calculatePotential(MultiVectorRCP conductance,
                                       VectorRCP sourceTerm) {

    MultiVectorRCP coordsMv = _coords->multiVector();
    MultiVectorRCP coefficients = _calculateCoefficients(conductance, coordsMv);
    MatrixRCP A = _buildGrid(coordsMv, coefficients);

    using precType = Ifpack2::Preconditioner<double, int, long long>;

    Ifpack2::Factory factory;
    Teuchos::RCP<precType> prec =
        factory.create<Tpetra::CrsMatrix<double, int, long long>>("ILUT", A);

    // TODO: Figure out the preconditioner parameters
    Teuchos::ParameterList precParams;
    precParams.set("fact: ilut level-of-fill", 2.0);
    precParams.set("fact: drop tolerance", 1e-4);
    prec->setParameters(precParams);
    prec->initialize();
    prec->compute();
    auto x = Teuchos::rcp(new Tpetra::Vector<double, int, long long>(_map));
    x->putScalar(0.0);

    auto rhs = Teuchos::rcp(
        new Tpetra::Vector<double, int, long long>(*sourceTerm, Teuchos::Copy));
    rhs->scale(RADIUS_EARTH_2);

    // TODO: Testing the pin point right now
    long long pinPoint = _nTh / 2;
    if (rhs->getMap()->isNodeGlobalElement(pinPoint)) {
        auto rhsData = rhs->getDataNonConst();
        auto localIdx = rhs->getMap()->getLocalElement(pinPoint);
        rhsData[localIdx] = 0.0;
    }

    auto problem =
        Teuchos::rcp(new Belos::LinearProblem<
                     double, Tpetra::MultiVector<double, int, long long>,
                     Tpetra::Operator<double, int, long long>>(A, x, rhs));

    problem->setRightPrec(prec);

    if (!problem->setProblem()) {
        throw std::runtime_error("Failed to set up problem");
    }

    auto solverParams = Teuchos::parameterList();
    solverParams->set("Maximum Iterations", 5000);
    solverParams->set("Convergence Tolerance", 1e-6);
    solverParams->set("Estimate Condition Number", true);

    int verbosity = Belos::Errors + Belos::Warnings + Belos::IterationDetails +
                    Belos::FinalSummary + Belos::TimingDetails +
                    Belos::StatusTestDetails;

    solverParams->set("Verbosity", verbosity);
    solverParams->set("Output Style", Belos::Brief);
    solverParams->set("Output Frequency", 10);

    Belos::BlockGmresSolMgr<double, Tpetra::MultiVector<double, int, long long>,
                            Tpetra::Operator<double, int, long long>>
        solver(problem, solverParams);

    Belos::ReturnType result = solver.solve();

    if (result == Belos::Converged) {
    } else {
        std::cout << "Failed to converge" << std::endl;
    }

    return x;
}

MatrixRCP TlSolver::_buildGrid(MultiVectorRCP coords,
                               MultiVectorRCP coefficients) {

    using std::sin;

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

        // TODO: Add a gauge condition at nTh / 2 (equator)

        // INPR: Replace with polar cap stencil
        // Applying the polar cap flux boundary condition
        if (theta == 0) {
            if (phi == 0) {
                // Our actual pole point
                matrixIdcs.push_back(gridPoint);
                vals.push_back(-sin(coords->dTh));
            } else {
                // Just constrain to phi == 0
                matrixIdcs.push_back(0);
                vals.push_back(1);
            }

        } else if (theta == _nTh - 1) {
        } else {
            matrixIdcs.push_back(gridPoint);
            vals.push_back(-2 * ththCoefficients[i] / (_dTh * _dTh) -
                           2 * phphCoefficients[i] / (_dPh * _dPh));

            size_t left = gridPoint - _nTh;
            size_t right = gridPoint + _nTh;

            size_t up = gridPoint - 1;
            size_t down = gridPoint + 1;

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
    }
    A->fillComplete();

    return A;
}

MultiVectorRCP TlSolver::_calculateCoefficients(MultiVectorRCP conductance,
                                                MultiVectorRCP coords) {
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
            up = oppositePhi * _nTh + 1;
        } else if (theta == _nTh - 1) {
            size_t oppositePhi = (phi + _nPh / 2) % _nPh;
            down = oppositePhi * _nTh + theta - 1;
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
