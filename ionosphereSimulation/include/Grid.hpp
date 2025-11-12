#pragma once
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <vector>

enum class Coords { TH, PH, COUNT };
enum class Ang { TH, PH, COUNT };
enum class Cart { X, Y, Z, COUNT };
enum class Sigma { THTH, THPH, PHPH, COUNT };
enum class DSigma { DTHTH_TH, DTHPH_PH, DTHPH_TH, DPHPH_PH, COUNT };
enum class HppSigma { HALL, PEDERSON, PARALLEL, COUNT };

template <typename E>
concept ValidEnum = std::is_same_v<E, Coords> || std::is_same_v<E, Sigma> ||
                    std::is_same_v<E, DSigma> || std::is_same_v<E, HppSigma> ||
                    std::is_same_v<E, Ang> || std::is_same_v<E, Cart>;

class Grid {
  public:
    Grid(size_t nTh, size_t nPh, double initialValue);
    double& operator()(size_t th, size_t ph);
    void printGrid();
    friend std::ostream& operator<<(std::ostream& outs, const Grid& g) {
        for (size_t th = 0; th < g.nTh; th++) {
            for (size_t ph = 0; ph < g.nPh; ph++) {
                outs << g._grid[th * g.nPh + ph];
                outs << ", ";
            }
            outs << "\n";
        }

        return outs;
    }

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

    friend std::ostream& operator<<(std::ostream& outs, const GridSet<T>& gs) {
        for (int grid = 0; grid < static_cast<int>(T::COUNT); grid++) {
            T test = static_cast<T>(grid);
            outs << *(gs._grids.at(test));
            outs << "\n";
        }
        return outs;
    }

    const size_t nTh;
    const size_t nPh;

  private:
    // TODO: Does this need to be a unique ptr
    std::map<T, std::unique_ptr<Grid>> _grids;
};
