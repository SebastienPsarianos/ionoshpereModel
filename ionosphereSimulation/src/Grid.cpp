#include "Grid.hpp"
#include <stdexcept>

Grid::Grid(size_t nTh, size_t nPh) : nTh(nTh), nPh(nPh), _grid(nTh * nPh, 0) {}

double& Grid::operator()(size_t th, size_t ph) {
    if (th >= nTh || ph >= nPh) {
        throw std::out_of_range("Index out of range");
    }
    return _grid[th * nTh + ph];
}

GridSet::GridSet(size_t count, size_t nTh, size_t nPh)
    : nTh(nTh), nPh(nPh), count(count) {
    for (size_t i = 0; i < count; i++) {
        _grids.push_back(std::make_unique<Grid>(nTh, nPh));
    }
}

double& GridSet::operator()(unsigned int idx, unsigned int th,
                            unsigned int ph) {
    if (th >= nTh) {
        throw std::out_of_range("Index THETA out of range");
    }
    if (ph >= nPh) {
        throw std::out_of_range("Index PHI out of range");
    }
    if (idx >= count) {
        throw new std::out_of_range("Grid index out of range");
    }
    return (*_grids[idx])(th, ph);
}
