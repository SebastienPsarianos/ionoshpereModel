#pragma once
#include <map>
#include <memory>
#include <optional>
#include <vector>

enum class Coords { TH, PH, COUNT };
enum class Sigma { THTH, THPH, PHPH, COUNT };
enum class DSigma { DTHTH_TH, DTHPH_PH, DTHPH_TH, DPHPH_PH, COUNT };
enum class HppSigma { HALL, PEDERSON, PARALLEL, COUNT };

template <typename E>
concept ValidEnum = std::is_same_v<E, Coords> || std::is_same_v<E, Sigma> ||
                    std::is_same_v<E, DSigma> || std::is_same_v<E, HppSigma>;

class Grid {
  public:
    Grid(size_t nTh, size_t nPh, double initialValue);
    double& operator()(size_t th, size_t ph);

    const size_t nTh;
    const size_t nPh;

  private:
    std::vector<double> _grid;
};

template <ValidEnum T> class GridSet {
  public:
    GridSet(size_t nTh, size_t nPh,
            std::optional<std::map<T, double>> initialValues = std::nullopt);
    double& operator()(T idx, unsigned int th, unsigned int ph);

    const size_t nTh;
    const size_t nPh;

  private:
    // TODO: Does this need to be a unique ptr
    std::map<T, std::unique_ptr<Grid>> _grids;
};
