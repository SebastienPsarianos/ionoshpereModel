#pragma once
#include <iostream>
#include <map>
#include <optional>
#include <vector>

enum class Coords { TH, PH, COUNT };
enum class Ang { TH, PH, COUNT };
enum class Cart { X, Y, Z, COUNT };
enum class Sigma { THTH, THPH, PHPH, COUNT };
enum class DSigma { DTHTH_TH, DTHPH_PH, DTHPH_TH, DPHPH_PH, COUNT };
enum class HppSigma { HALL, PEDERSON, PARALLEL, COUNT };
enum class Coeff { THTH, PHPH, TH, PH, COUNT };

struct {
    double value;
    double th_val;
    double ph_val;
} typedef GridVal;

template <typename E>
concept ValidEnum = std::is_same_v<E, Coords> || std::is_same_v<E, Sigma> ||
                    std::is_same_v<E, DSigma> || std::is_same_v<E, HppSigma> ||
                    std::is_same_v<E, Ang> || std::is_same_v<E, Cart> ||
                    std::is_same_v<E, Coeff>;

template <ValidEnum T> class GridSet;

class Grid {
  public:
    Grid(size_t nTh, size_t nPh, double initialValue);
    double& operator()(size_t th, size_t ph);
    std::ostream& printWithCoords(std::ostream& out, GridSet<Ang>& coords);

    friend std::ostream& operator<<(std::ostream& outs, const Grid& g) {
        for (size_t th = 0; th < g.nTh; th++) {
            for (size_t ph = 0; ph < g.nPh; ph++) {
                outs << g._grid.at(th * g.nPh + ph);
                outs << ", ";
            }
            outs << "\n";
        }

        return outs;
    }

    size_t nTh;
    size_t nPh;

  private:
    std::vector<double> _grid;
};

template <ValidEnum T> class GridSet {
  public:
    GridSet(size_t nTh, size_t nPh,
            std::optional<std::map<T, double>> initialValues = std::nullopt);
    double& operator()(T idx, unsigned int th, unsigned int ph);
    std::ostream& printWithCoords(std::ostream& out, GridSet<Ang>& coords);

    friend std::ostream& operator<<(std::ostream& outs, const GridSet<T>& gs) {
        for (int grid = 0; grid < static_cast<int>(T::COUNT); grid++) {
            T test = static_cast<T>(grid);
            outs << gs._grids.at(test);
            outs << "\n";
        }
        return outs;
    }

    size_t nTh;
    size_t nPh;

  private:
    std::map<T, Grid> _grids;
};
