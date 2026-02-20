#include "ionosphere/solver/Solver.hpp"
#include "ionosphere/utils/Grid.hpp"

class TlSolver : Solver {
  public:
    TlSolver(size_t nTh, size_t nPh, double dPh, double dTh, MapRcp map);

    // TODO: Fix solver
    void calculatePotential(MultiVectorRcp radCurrent,
                            MultiVectorRcp conductance, MultiVectorRcp coords);

  private:
    void _buildGrid(MultiVectorRcp radCurrent, MultiVectorRcp conductance,
                    MultiVectorRcp coords);

    void _calculateCoefficients(MultiVectorRcp conductance,
                                MultiVectorRcp coords);

    size_t _nTh;
    size_t _nPh;
    double _dTh;
    double _dPh;

    MapRcp _map;
    MultiVectorRcp _coefficients;
};
