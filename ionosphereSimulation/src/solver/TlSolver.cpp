#include "ionosphere/solver/TlSolver.hpp"
#include "ionosphere/TrilinosAliases.hpp"

#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Operator.hpp>
#include <iostream>
#include <stdexcept>

Ionosphere::VectorRCP
Ionosphere::calculatePotential(Teuchos::RCP<Problem> problem) {
    using precType = Ifpack2::Preconditioner<Scalar, LocalOrd, GlobalOrd>;

    using namespace Ionosphere;

    // Build the problem that was passed in
    problem->build();

    Ifpack2::Factory factory;
    Teuchos::RCP<precType> prec =
        factory.create<Matrix>("ILUT", problem->matrix);

    Teuchos::ParameterList precParams;
    precParams.set("fact: ilut level-of-fill", 2.0);
    precParams.set("fact: drop tolerance", 1e-4);
    prec->setParameters(precParams);
    prec->initialize();
    prec->compute();

    auto belosProblem = Teuchos::rcp(
        new Belos::LinearProblem<Scalar, MultiVector,
                                 Tpetra::Operator<Scalar, LocalOrd, GlobalOrd>>(
            problem->matrix, problem->guess, problem->rhs));

    belosProblem->setRightPrec(prec);

    if (!belosProblem->setProblem()) {
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
        solver(belosProblem, solverParams);

    Belos::ReturnType result = solver.solve();

    if (result == Belos::Converged) {
    } else {
        std::cout << "Failed to converge" << std::endl;
    }

    return problem->guess;
}
