#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/conductance/EmpConductance.hpp"
#include "ionosphere/coordinates/Coordinates.hpp"
#include "ionosphere/io/tecplot.hpp"
#include "ionosphere/mhd/MHDInterpolation.hpp"
#include "ionosphere/solver/Problem.hpp"
#include "ionosphere/solver/ProblemFactory.hpp"
#include "ionosphere/solver/TlSolver.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_YAML.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Tpetra_Core.hpp>
#include <stdexcept>

int main(int argc, char* argv[]) {

    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    auto comm = Tpetra::getDefaultComm();

    /** START PARAMETER PARSING **/

    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromYamlFileAndBroadcast(
        "configuration.yaml", Teuchos::Ptr<Teuchos::ParameterList>(&paramList),
        *comm);

    auto conductanceParams = paramList.sublist("conductance");
    auto ioParams = paramList.sublist("io");
    auto equationParams = paramList.sublist("equation");
    auto solverParams = paramList.sublist("solver");

    /*** END PARAMETER PARSING ***/

    /*** BEGIN MHD INTERPOLATION ***/
    int nTh = -1;
    int nPh = -1;

    MHDInterpolation::getGridSize(ioParams, &nTh, &nPh, comm);
    if (nTh <= 0 || nPh <= 0) {
        throw std::runtime_error(
            "Error reading source term data, ending execution");
    }

    auto map = Teuchos::rcp(new Ionosphere::Map(nTh * nPh, 0, comm));

    auto sourceTerm = Teuchos::rcp(new Ionosphere::Vector(map));

    Teuchos::RCP<Coordinates> coords;

    auto solarModel =
        Teuchos::rcp<SolarModel>(new SolarModel(conductanceParams));

    auto dipoleModel =
        Teuchos::rcp<DipoleModel>(new DipoleModel(comm, conductanceParams));

    MHDInterpolation::processLegacyOutput(coords, dipoleModel, solarModel,
                                          sourceTerm, map, comm, nTh, nPh,
                                          ioParams);
    /*** END MHD INTERPOLATION ***/

    /*** BEGIN CONDUCTANCE CALC ***/
    auto [auroralConductance, euvConductance, conductance] =
        EmpConductance(coords, map, conductanceParams).computeConductance();
    /*** END CONDUCTANCE CALC ***/

    /*** BEGIN POTENTIAL SOLVE ***/
    Teuchos::RCP<Problem> problem = Ionosphere::problemFactory(
        equationParams, coords, sourceTerm, conductance, map);

    Ionosphere::VectorRCP result = Ionosphere::calculatePotential(problem);
    /*** END POTENTIAL SOLVE ***/

    /*** BEGIN PLOTTING ***/
    Ionosphere::exportToTecplot(ioParams, coords, result, sourceTerm,
                                auroralConductance, euvConductance, conductance,
                                comm, coords->nTh, coords->nPh);
    /*** END PLOTTING ***/

    return 0;

    // TODO: Proper post-processing, including J and others
}
