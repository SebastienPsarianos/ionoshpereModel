#include "Grid.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

Grid::Grid(size_t nTh, size_t nPh, double initialValue)
    : nTh(nTh), nPh(nPh), _grid(nTh * nPh, initialValue) {}

double& Grid::operator()(size_t th, size_t ph) {
    if (th >= nTh || ph >= nPh) {
        throw std::out_of_range("Index out of range");
    }
    return _grid[th * nPh + ph];
}

// TODO: Redo this properly.
std::ostream& Grid::printWithCoords(std::ostream& out, GridSet<Ang>& coords) {
    out << std::scientific << std::setprecision(6);

    for (size_t th = 0; th < nTh; th++) {
        for (size_t ph = 0; ph < nPh; ph++) {
            out << coords(Ang::TH, th, ph) << " " << coords(Ang::PH, th, ph)
                << " " << (*this)(th, ph) << "\n";
        }
    }
    out << std::endl;
    return out;
}

template <ValidEnum T>
GridSet<T>::GridSet(size_t nTh, size_t nPh,
                    std::optional<std::map<T, double>> initialValues)
    : nTh(nTh), nPh(nPh) {

    if (!initialValues.has_value()) {
        for (int grid = 0; grid < static_cast<int>(T::COUNT); grid++) {
            _grids.try_emplace(static_cast<T>(grid), nTh, nPh, 0);
        }
        return;
    }

    if (initialValues.value().size() != (int)T::COUNT) {
        std::cerr << "Initial values not defined";
        // TODO: Pick a better exception
        throw std::out_of_range("Initial values not defined for all entries");
    }

    for (const auto& [grid, value] : initialValues.value()) {
        _grids.try_emplace(static_cast<T>(grid), nTh, nPh, value);
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
    return _grids.at(grid)(th, ph);
}

template <ValidEnum T>
std::ostream& GridSet<T>::printWithCoords(std::ostream& out,
                                          GridSet<Ang>& coords) {
    out << std::scientific << std::setprecision(6);

    for (size_t th = 0; th < nTh; th++) {
        for (size_t ph = 0; ph < nPh; ph++) {
            out << coords(Ang::TH, th, ph) << " " << coords(Ang::PH, th, ph)
                << " ";

            // Print out values for all of the grid values in a row
            for (int grid = 0; grid < static_cast<int>(T::COUNT); grid++) {
                out << this->_grids.at(static_cast<T>(grid))(th, ph) << " ";
            }

            out << "\n";
        }
    }

    out << std::endl;
    return out;
}

template class GridSet<Sigma>;
template class GridSet<Ang>;
template class GridSet<Coords>;
template class GridSet<DSigma>;
template class GridSet<HppSigma>;
template class GridSet<Cart>;
template class GridSet<Coeff>;
