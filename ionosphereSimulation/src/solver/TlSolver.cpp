#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/TrilinosAliases.hpp"

#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Operator.hpp>
#include <cmath>
#include <stdexcept>

// TODO: Come up with a better setup
const long long IONOSPHERE_Radius_Earth = 6378000.00;
const long long IONOSPHERE_Height_Earth = 400000.00;
const long long RADIUS_EARTH =
    (IONOSPHERE_Height_Earth + IONOSPHERE_Radius_Earth);
const long long RADIUS_EARTH_2 = (RADIUS_EARTH * RADIUS_EARTH);

using namespace Ionosphere;

TlSolver::TlSolver(Teuchos::RCP<Coordinates> coords,
                   Ionosphere::MultiVectorRCP conductance,
                   Ionosphere::VectorRCP sourceTerm, Ionosphere::MapRCP map)
    : _coords(coords), _conductance(conductance), _sourceTerm(sourceTerm),
      _map(map) {}

VectorRCP TlSolver::calculatePotential() {

    MultiVectorRCP coefficients = _calculateCoefficients();
    MatrixRCP A = _buildGrid(coefficients);

    using precType = Ifpack2::Preconditioner<Scalar, LocalOrd, GlobalOrd>;

    Ifpack2::Factory factory;
    Teuchos::RCP<precType> prec = factory.create<Matrix>("ILUT", A);

    // TODO: Do some more preconditioner and solver optimization
    //  - Look into zoltan for map generation
    Teuchos::ParameterList precParams;
    precParams.set("fact: ilut level-of-fill", 2.0);
    precParams.set("fact: drop tolerance", 1e-4);
    prec->setParameters(precParams);
    prec->initialize();
    prec->compute();
    auto x = Teuchos::rcp(new Vector(_map));
    x->putScalar(0.0);

    auto rhs = Teuchos::rcp(new Vector(*_sourceTerm, Teuchos::Copy));
    rhs->scale(RADIUS_EARTH_2);

    GlobalOrd pinPoint = _coords->nTh / 2;
    {
        auto rhsData = rhs->getDataNonConst();
        auto myGridPoints = _map->getMyGlobalIndices();
        for (LocalOrd i = 0; i < static_cast<LocalOrd>(myGridPoints.size());
             i++) {
            GlobalOrd currentGid = myGridPoints[i];
            GlobalOrd theta = currentGid % _coords->nTh;
            GlobalOrd phi = currentGid / _coords->nTh;

            if (currentGid == pinPoint) {
                rhsData[i] = 0.0;
            } else if (theta == 0 || theta == _coords->nTh - 1) {
                // We set the two primary pole points to the cap condition
                if (phi == 0) {
                    Scalar theta0 = _coords->dTh / 2.0;
                    Scalar jR_pole = rhsData[i] / RADIUS_EARTH_2;
                    rhsData[i] = jR_pole * RADIUS_EARTH_2 * 2.0 * M_PI *
                                 (1.0 - std::cos(theta0));
                } else {
                    // we constrain the rest of the pole points to theta
                    rhsData[i] = 0.0;
                }
            }
        }
    }

    auto problem = Teuchos::rcp(
        new Belos::LinearProblem<Scalar, MultiVector,
                                 Tpetra::Operator<Scalar, LocalOrd, GlobalOrd>>(
            A, x, rhs));

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

    Belos::BlockGmresSolMgr<Scalar, MultiVector,
                            Tpetra::Operator<Scalar, LocalOrd, GlobalOrd>>
        solver(problem, solverParams);

    Belos::ReturnType result = solver.solve();

    if (result == Belos::Converged) {
    } else {
        std::cout << "Failed to converge" << std::endl;
    }

    return x;
}

MultiVectorRCP TlSolver::_gatherPoleData() {
    auto comm = _map->getComm();
    size_t localSize = 2 * (_coords->nPh + 1);

    // Create a local map that contains all pole values
    auto localMap = Teuchos::rcp(
        new Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                localSize, 0, comm));

    // 0 = theta, 1 = sigThTh, 2 = sigThPh (needed for the pole boundary
    // condition)
    auto poleData = Teuchos::rcp(new MultiVector(localMap, 3));
    poleData->putScalar(0.0);

    auto pdTh = poleData->getDataNonConst(0);
    auto sigThTh = poleData->getDataNonConst(1);
    auto sigThPh = poleData->getDataNonConst(2);

    auto coordsMv = _coords->multiVector();
    auto thVals = coordsMv->getDataNonConst(0);
    auto pedVals = _conductance->getDataNonConst(0);
    auto hallVals = _conductance->getDataNonConst(1);
    auto parVals = _conductance->getDataNonConst(2);

    // Construct lookup table of the entries we want to keep
    // 0 is the pole and then the rest are the ring nodes. pdIdx gives the index
    // in the localVector. TODO: This needs to be cleaner
    struct PoleEntry {
        GlobalOrd currentGid;
        GlobalOrd pdIdx;
    };

    std::vector<PoleEntry> entries;

    // Add the north pole
    entries.push_back({0, 0});
    // Add all of the north pole ring points
    for (GlobalOrd j = 0; j < _coords->nPh; j++)
        entries.push_back(
            {static_cast<GlobalOrd>(j * _coords->nTh + 1), j + 1});

    // Add the south pole
    entries.push_back({_coords->nTh - 1, _coords->nPh + 1});

    // Add the south pole ring points
    for (GlobalOrd j = 0; j < _coords->nPh; j++)
        entries.push_back(
            {j * _coords->nTh + _coords->nTh - 2, _coords->nPh + 2 + j});

    // Compute the required spherical conductance values for each pole point
    for (const auto& entry : entries) {
        if (!_map->isNodeGlobalElement(entry.currentGid))
            continue;

        auto li = _map->getLocalElement(entry.currentGid);
        Scalar th = thVals[li];
        Scalar cosE = (-2.0 * std::cos(th)) /
                      std::sqrt(1.0 + 3.0 * std::cos(th) * std::cos(th));
        Scalar sinE =
            std::sin(th) / std::sqrt(1.0 + 3.0 * std::cos(th) * std::cos(th));
        Scalar C = parVals[li] * cosE * cosE + pedVals[li] * sinE * sinE;

        pdTh[entry.pdIdx] = th;
        sigThTh[entry.pdIdx] = parVals[li] * pedVals[li] / C;
        sigThPh[entry.pdIdx] = parVals[li] * hallVals[li] * (-cosE) / C;
    }

    // Use a reduce to get the values available in all processors
    auto poleDataLocal = Teuchos::rcp(new MultiVector(*poleData));
    for (size_t col = 0; col < 3; col++) {
        auto src = poleDataLocal->getDataNonConst(col);
        auto dst = poleData->getDataNonConst(col);
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM,
                           static_cast<int>(localSize), src.getRawPtr(),
                           dst.getRawPtr());
    }
    return poleData;
}
MatrixRCP TlSolver::_buildGrid(MultiVectorRCP coefficients) {

    using std::sin;

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

    auto poleData = _gatherPoleData();
    auto sigThTh = poleData->getData(1);
    auto sigThPh = poleData->getData(2);

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

        GlobalOrd primaryPoleId = (theta == 0) ? 0 : (_coords->nTh - 1);

        // Applying the polar cap flux boundary condition (see notes)
        if (theta == 0 || theta == _coords->nTh - 1) {
            if (phi == 0) {
                Scalar theta0 = _coords->dTh / 2;

                // Grab the indecies in the poleData multivector
                LocalOrd pole = theta == 0 ? 0 : (_coords->nPh + 1);
                LocalOrd ringStart = (theta == 0) ? 1 : (_coords->nPh + 2);

                // Multiply terms by negative 1 at the south pole to
                // account for  negative normal vector direction
                Scalar multiplier = theta == 0 ? 1 : -1;

                // Summing average conductance contributions across all phi
                // for pole cap (taking an average of the pole and the first
                // ring of theta values). Here we are using the 0 point for
                // consistency sake.
                Scalar poleCoefficient = 0;
                for (GlobalOrd i = 0; i < _coords->nPh; i++) {
                    poleCoefficient +=
                        (sigThTh[ringStart + i] + sigThTh[pole]) / 2;
                }
                poleCoefficient *= -sin(theta0) * _coords->dPh / _coords->dTh;

                matrixIdcs.push_back(gridPoint);
                vals.push_back(poleCoefficient);

                // Handle the contributions from the surrounding ring

                for (GlobalOrd i = 0; i < _coords->nPh; i++) {
                    Scalar sigThThAvg =
                        (sigThTh[ringStart + i] + sigThTh[pole]) / 2;

                    GlobalOrd prevIdx = (i - 1 + _coords->nPh) % _coords->nPh;
                    GlobalOrd nextIdx = (i + 1) % _coords->nPh;
                    Scalar sigThPhAvg = (sigThPh[ringStart + prevIdx] -
                                         sigThPh[ringStart + nextIdx]) /
                                        8;

                    Scalar ringPointCoeff =

                        (_coords->dPh * sin(theta0) / _coords->dTh) *
                            sigThThAvg +
                        multiplier * sigThPhAvg;

                    GlobalOrd ringTheta = (theta == 0) ? 1 : (_coords->nTh - 2);
                    GlobalOrd ringGlobIdx = i * _coords->nTh + ringTheta;

                    matrixIdcs.push_back(ringGlobIdx);
                    vals.push_back(ringPointCoeff);
                }

                A->insertGlobalValues(gridPoint, matrixIdcs, vals);

            } else {
                // Just constrain to phi == 0
                matrixIdcs.push_back(gridPoint);
                vals.push_back(1);
                matrixIdcs.push_back(primaryPoleId);
                vals.push_back(-1);
                A->insertGlobalValues(gridPoint, matrixIdcs, vals);
            }
        } else {
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
    A->fillComplete();

    return A;
}

MultiVectorRCP TlSolver::_calculateCoefficients() {
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
