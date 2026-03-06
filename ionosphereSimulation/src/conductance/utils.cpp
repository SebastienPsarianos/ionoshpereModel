#include "ionosphere/conductance/utils.hpp"

#include <cmath>

double fourierSeries(nlohmann::json coefficients, double mlt) {
    using std::cos, std::sin;
    double coefficient = coefficients["const"].get<double>();
    for (int i = 0; i < 6; i++) {
        coefficient += coefficients["cos"][i].get<double>() *
                           cos((i + 1) * mlt * M_PI / 12) +
                       coefficients["sin"][i].get<double>() *
                           sin((i + 1) * mlt * M_PI / 12);
    }
    return coefficient;
}

double epsteinFunction(double h, double r, double h0, double S1, double S2) {
    return r + S1 * (h - h0) +
           (S2 - S1) * std::log((1.0 - (S1 / (S2 * std::exp(-(h - h0))))) /
                                (1.0 - (S1 / S2)));
}
