#include "ionosphere/utils/Grid.hpp"
#include <optional>
#include <stdexcept>

template <typename T>
Grid<T>::Grid(size_t nTh, size_t nPh, std::optional<T> initialValue)
    : nTh(nTh), nPh(nPh), _grid(nTh * nPh, initialValue.value_or(T())) {}

template <typename T> T& Grid<T>::operator()(size_t th, size_t ph) {
    if (th >= nTh || ph >= nPh) {
        throw std::out_of_range("Index out of range");
    }
    return _grid[th * nPh + ph];
}
template <typename T> const T& Grid<T>::operator()(size_t th, size_t ph) const {
    if (th >= nTh || ph >= nPh) {
        throw std::out_of_range("Index out of range");
    }
    return _grid[th * nPh + ph];
}

template <typename T>
std::ostream& Grid<T>::printWithCoords(std::ostream& out,
                                       const Grid<ThPh>& coords) {
    for (size_t th = 0; th < nTh; th++) {
        for (size_t ph = 0; ph < nPh; ph++) {
            out << coords(th, ph).th << " " << coords(th, ph).ph << " "
                << (*this)(th, ph) << "\n";
        }
    }
    return out;
}

template class Grid<double>;
template class Grid<ThPh>;
template class Grid<CartVector>;
template class Grid<Sigma>;
template class Grid<HppSigma>;
template class Grid<DSigma>;
template class Grid<Coeff>;
