#include "Grid.hpp"
#include <stdexcept>

Grid::Grid(size_t nTh, size_t nPh, double initialValue)
    : nTh(nTh), nPh(nPh), _grid(nTh * nPh, initialValue) {}

double& Grid::operator()(size_t th, size_t ph) {
    if (th >= nTh || ph >= nPh) {
        throw std::out_of_range("Index out of range");
    }
    return _grid[th * nPh + ph];
}

template <ValidEnum T>
GridSet<T>::GridSet(size_t nTh, size_t nPh,
                    std::optional<std::map<T, double>> initialValues)
    : nTh(nTh), nPh(nPh) {

    if (!initialValues.has_value()) {
        for (int i = 0; i < static_cast<int>(T::COUNT); i++) {
            _grids[static_cast<T>(i)] = std::make_unique<Grid>(nTh, nPh, 0);
        }
        return;
    }

    if (initialValues.value().size() != (int)T::COUNT) {
        std::cerr << "Initial values not defined";
        // TODO: Pick a better exception
        throw std::out_of_range("Initial values not defined for all entries");
    }

    for (const auto& [grid, value] : initialValues.value()) {
        _grids[grid] = std::make_unique<Grid>(nTh, nPh, value);
    }
}

template <ValidEnum T>
double& GridSet<T>::operator()(T grid, unsigned int th, unsigned int ph) {
    if (th >= nTh) {
        throw std::out_of_range("Index THETA out of range");
    }
    if (ph >= nPh) {
        throw std::out_of_range("Index PHI out of range");
    }
    return (*_grids[grid])(th, ph);
}

template class GridSet<Sigma>;
template class GridSet<Ang>;
template class GridSet<Coords>;
template class GridSet<DSigma>;
template class GridSet<HppSigma>;
template class GridSet<Cart>;
