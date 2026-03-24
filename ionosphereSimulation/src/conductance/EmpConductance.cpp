#include "ionosphere/conductance/EmpConductance.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ionosphere/TrilinosAliases.hpp"
#include "ionosphere/conductance/EuvConductance.hpp"
#include "ionosphere/conductance/HardyConductance.hpp"

#include <cmath>
#include <nlohmann/json>

using namespace Ionosphere;
using nlohmann::json;

EmpConductance::EmpConductance(Teuchos::RCP<Coordinates> coords, MapRCP map,
                               const Teuchos::ParameterList& conductanceParams)
    : _map(map), _coords(coords), hardyConductanceModel(coords, map),
      euvConductanceModel(coords, map), _kp(conductanceParams.get<int>("kp")),
      _f107(conductanceParams.get<double>("f107")),
      _sig0(conductanceParams.get<double>("sig0")) {}

std::tuple<MultiVectorRCP, MultiVectorRCP, MultiVectorRCP>
EmpConductance::computeConductance() {
    MultiVectorRCP auroralConductance =
        hardyConductanceModel.computeAuroralConductance(_kp);
    MultiVectorRCP euvConductance =
        euvConductanceModel.computeEuvConductance(_f107);
    MultiVectorRCP conductance =
        _computeHppConductance(auroralConductance, euvConductance);

    return {auroralConductance, euvConductance, conductance};
}

MultiVectorRCP
EmpConductance::_computeHppConductance(MultiVectorRCP auroralConductance,
                                       MultiVectorRCP euvConductance) {
    MultiVectorRCP conductance =
        Teuchos::rcp(new Ionosphere::MultiVector(_map, 3));

    auto pedersonConductancesAur = auroralConductance->getDataNonConst(0);
    auto hallConductancesAur = auroralConductance->getDataNonConst(1);
    auto pedersonConductancesEuv = euvConductance->getDataNonConst(0);
    auto hallConductancesEuv = euvConductance->getDataNonConst(1);

    auto hppPedersonConductance = conductance->getDataNonConst(0);
    auto hppHallConductance = conductance->getDataNonConst(1);
    auto hppParallelConductance = conductance->getDataNonConst(2);

    for (size_t i = 0; i < _coords->multiVector()->getLocalLength(); i++) {
        hppHallConductance[i] =
            sqrt(hallConductancesAur[i] * hallConductancesAur[i] +
                 hallConductancesEuv[i] * hallConductancesEuv[i]);
        hppPedersonConductance[i] =
            sqrt(pedersonConductancesAur[i] * pedersonConductancesAur[i] +
                 pedersonConductancesEuv[i] * pedersonConductancesEuv[i]);

        // Imposing a minimum background conductance of 0.25 mhos
        hppPedersonConductance[i] = std::max(hppPedersonConductance[i], 0.25);
        hppHallConductance[i] = std::max(hppHallConductance[i], 0.25);
        hppParallelConductance[i] = _sig0;
    }
    return conductance;
}
