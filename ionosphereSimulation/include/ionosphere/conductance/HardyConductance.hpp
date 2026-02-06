#pragma once
#include <array>
#include <cmath>

struct FourierCoeffs {
    double constant;
    std::array<double, 6> cos_terms;
    std::array<double, 6> sin_terms;
};

struct EpsteinParams {
    FourierCoeffs max_val;    // r
    FourierCoeffs max_lat;    // h0
    FourierCoeffs up_slope;   // s1
    FourierCoeffs down_slope; // s2
};

// This tells the compiler: "This array exists in HardyCoefficients.cpp"
extern const std::array<EpsteinParams, 7> HARDY_PEDERSON_COEFFS;
extern const std::array<EpsteinParams, 7> HARDY_HALL_COEFFS;
