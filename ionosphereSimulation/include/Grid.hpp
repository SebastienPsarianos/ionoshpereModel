#pragma once
#include <vector>

class Grid {
  public:
    Grid(size_t nTh, size_t nPh);
    double& operator()(size_t th, size_t ph);
    const size_t nTh;
    const size_t nPh;

  private:
    std::vector<double> _grid;
};

class GridSet {
  public:
    GridSet(size_t count, size_t nTh, size_t nPh);
    double& operator()(unsigned int idx, unsigned int th, unsigned int ph);

    const size_t nTh;
    const size_t nPh;
    const size_t count;

  private:
    std::vector<std::unique_ptr<Grid>> _grids;
};
