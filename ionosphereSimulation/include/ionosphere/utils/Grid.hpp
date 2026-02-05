#pragma once
#include <iostream>
#include <optional>
#include <vector>

struct {
    double th;
    double ph;
} typedef ThPh;

struct {
    double x;
    double y;
    double z;
} typedef CartVector;

struct {
    double thth;
    double thph;
    double phph;
} typedef Sigma;

struct {
    double hall;
    double pederson;
    double parallel;
} typedef HppSigma;

struct {
    double dthth_th;
    double dthph_ph;
    double dthph_th;
    double dphph_ph;
} typedef DSigma;

struct {
    double thth;
    double phph;
    double th;
    double ph;
} typedef Coeff;

template <typename T>
concept ValidStruct =
    std::is_same_v<T, ThPh> || std::is_same_v<T, CartVector> ||
    std::is_same_v<T, Sigma> || std::is_same_v<T, HppSigma> ||
    std::is_same_v<T, DSigma> || std::is_same_v<T, Coeff> ||
    std::is_same_v<T, double>;

template <typename T> class Grid {
  public:
    Grid(size_t nTh, size_t nPh, std::optional<T> initialValue = std::nullopt);
    T& operator()(size_t th, size_t ph);
    const T& operator()(size_t th, size_t ph) const;
    std::ostream& printWithCoords(std::ostream& out, const Grid<ThPh>& coords);

    size_t nTh;
    size_t nPh;

  private:
    std::vector<T> _grid;
};

inline std::ostream& operator<<(std::ostream& out, const ThPh& s) {
    return out << s.th << " " << s.ph;
}

inline std::ostream& operator<<(std::ostream& out, const CartVector& s) {
    return out << s.x << " " << s.y << " " << s.z;
}

inline std::ostream& operator<<(std::ostream& out, const Sigma& s) {
    return out << s.thth << " " << s.thph << " " << s.phph;
}

inline std::ostream& operator<<(std::ostream& out, const HppSigma& s) {
    return out << s.hall << " " << s.pederson << " " << s.parallel;
}

inline std::ostream& operator<<(std::ostream& out, const DSigma& s) {
    return out << s.dthth_th << " " << s.dthph_ph << " " << s.dthph_th << " "
               << s.dphph_ph;
}

inline std::ostream& operator<<(std::ostream& out, const Coeff& s) {
    return out << s.thth << " " << s.phph << " " << s.th << " " << s.ph;
}
