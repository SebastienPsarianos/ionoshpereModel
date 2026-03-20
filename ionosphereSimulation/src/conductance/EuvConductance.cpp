#include "ionosphere/conductance/EuvConductance.hpp"
#include "ionosphere/TrilinosAliases.hpp"

#include <cmath>

using namespace Ionosphere;

EuvConductance::EuvConductance(Teuchos::RCP<Coordinates> coords, MapRCP map)
    : _map(map), _euvConductance(new Ionosphere::MultiVector(map, 3)),
      _coords(coords) {}

MultiVectorRCP EuvConductance::computeEuvConductance(double f107) {
    using std::cos, std::sin, std::sqrt;
    auto pedersonConductances = _euvConductance->getDataNonConst(0);
    auto hallConductances = _euvConductance->getDataNonConst(1);

    for (size_t i = 0; i < _coords->multiVector()->getLocalLength(); i++) {
        double sza = _coords->localIdx2Sza(i);

        if (sza >= M_PI / 2) {
            pedersonConductances[i] = 0;
            hallConductances[i] = 0;
        } else {
            pedersonConductances[i] =
                pow(f107, 0.49) * (0.34 * cos(sza) + 0.93 * sqrt(cos(sza)));
            hallConductances[i] =
                pow(f107, 0.53) * (0.81 * cos(sza) + 0.54 * sqrt(cos(sza)));
        }
    }

    return _euvConductance;
}
