#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/solver/Solver.hpp"
#include "ionosphere/utils/Coordinates.hpp"

class TlSolver : Solver {
  public:
    TlSolver(Teuchos::RCP<Coordinates> coords, Ionosphere::MapRCP map);

    Ionosphere::VectorRCP
    calculatePotential(Ionosphere::MultiVectorRCP conductance,
                       Ionosphere::VectorRCP sourceTerm);

  private:
    Ionosphere::MultiVectorRCP
    _calculateCoefficients(Ionosphere::MultiVectorRCP conductance,
                           Ionosphere::MultiVectorRCP coords);

    Ionosphere::VectorRCP
    _buildSourceVector(Ionosphere::MultiVectorRCP sourceTerm);

    Ionosphere::MatrixRCP _buildGrid(Ionosphere::MultiVectorRCP coords,
                                     Ionosphere::MultiVectorRCP coefficients);

    double _dTh;
    double _dPh;

    Teuchos::RCP<Coordinates> _coords;
    Ionosphere::MapRCP _map;
};
