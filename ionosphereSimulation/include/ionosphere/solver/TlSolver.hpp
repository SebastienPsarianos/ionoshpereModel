#include "ionosphere/IonosphereTypes.hpp"
#include "ionosphere/solver/Solver.hpp"

class TlSolver : Solver {
  public:
    TlSolver(size_t nTh, size_t nPh, double dPh, double dTh,
             Ionosphere::MapRCP map);

    // TODO: Fix solver
    Ionosphere::VectorRCP
    calculatePotential(Ionosphere::MultiVectorRCP conductance,
                       Ionosphere::MultiVectorRCP coords,
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

    Ionosphere::MapRCP _map;
};
