#include "Grid.hpp"

#define MAX_ITERATION_NUM 100000
#define RES_THRESHOLD 10

#define IONOSPHERE_Radius_Earth 6378000.00
#define IONOSPHERE_Height_Earth 400000.00
#define RADIUS_EARTH IONOSPHERE_Height_Earth + IONOSPHERE_Radius_Earth

enum Algorithm { GAUSS_SEIDEL, JACOBI, SOR };

class Solver {
  public:
    Solver(int nTh, int nPh, std::shared_ptr<GridSet<Coeff>> kappa,
           std::shared_ptr<GridSet<Ang>> coords,
           std::shared_ptr<Grid> radCurrent, Algorithm algorithm);
    std::shared_ptr<Grid> calculatePotential();

  private:
    size_t nTh;
    size_t nPh;
    double dTh;
    double dPh;
    double dTh2;
    double dPh2;
    std::shared_ptr<Grid> radCurrent;
    std::shared_ptr<GridSet<Ang>> coords;
    std::shared_ptr<GridSet<Coeff>> kappa;
    std::shared_ptr<Grid> potential;
    Algorithm algorithm;
    Grid previousIteration;

    double _gaussSeidelFormula(size_t th, size_t ph);
    double _calculateResidual(size_t th, size_t ph);
};
