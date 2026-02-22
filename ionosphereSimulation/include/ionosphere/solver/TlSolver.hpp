#include "ionosphere/solver/Solver.hpp"
#include "ionosphere/utils/Grid.hpp"

class TlSolver : Solver {
  public:
    TlSolver(size_t nTh, size_t nPh, double dPh, double dTh, MapRcp map);

    // TODO: Fix solver
    VectorRcp calculatePotential(MultiVectorRcp conductance,
                                 MultiVectorRcp coords, VectorRcp sourceTerm);

  private:
    MultiVectorRcp _calculateCoefficients(MultiVectorRcp conductance,
                                          MultiVectorRcp coords);

    VectorRcp _buildSourceVector(MultiVectorRcp sourceTerm);

    CrsMatrixRcp _buildGrid(MultiVectorRcp coords, MultiVectorRcp coefficients);

    double _dTh;
    double _dPh;

    MapRcp _map;
};
