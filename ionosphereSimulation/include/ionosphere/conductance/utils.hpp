#pragma once
#include <nlohmann/json>

double fourierSeries(nlohmann::json coefficients, double mlt);
double epsteinFunction(double h, double r, double h0, double S1, double S2);
