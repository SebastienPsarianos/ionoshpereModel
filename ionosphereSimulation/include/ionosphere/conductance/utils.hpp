#pragma once
#include "ionosphere/utils/Grid.hpp"

#include <Eigen/Dense>
#include <nlohmann/json>
#define SUNFREQ .0172019715

// TODO: Check over this whole thing and likely turn into class containing
// rotation matrix and  the time stuff

/**
 * @brief Calculates the solar zenith angle on the earth's surface using
 * algorithm 3 from Grena (2012)
 * TODO: rewrite docs
 *
 *
 * It applies a Parallax Correction to adjust for the observer's position on
 * the Earth's surface (Topocentric), but doesn't apply the refraction
 * correction.
 *
 * @param[in] utDays     UT time in days (see computeGrenaTimeScales).
 * @param[in] ttDays     UT time in days (see computeGrenaTimeScales).
 * @param[in] latitude   Latitude in radians (-PI/2 to +PI/2).
 * @param[in] longitude  Longitude in radians (East positive).
 *
 * @return The Solar Zenith Angle in radians.
 */
double computeSolarZenith(double utTime, double ttTime, GeoGeo coords);

/**
 * @brief Converts a standard date into the specific time scales
 * required by the Grena (2012) algorithms.
 *
 * This function calculates:
 * 1. UT Time (t): Days from the beginning of 2060 (Universal Time).
 * 2. TT Time (te): Days from the beginning of 2060 (Terrestrial Time)
 *
 * @note Valid for years 2010 to 2110.
 *
 * @param[out] utDays  (Output) Time scale based on UT. Corresponds to variable
 * 't' in Grena paper
 *  @param[out] ttDays  (Output) Time scale based on TT. Corresponds to variable
 * 'te' in Grena paper.
 *
 * @param[in]  year    Year 2010-2110).
 * @param[in]  month   Month (1-12).
 * @param[in]  day     Day of the month (1-31).
 * @param[in]  hour    Hour of the day (0-23).
 */
void computeGrenaTimescales(double& utTime, double& ttTime, int year, int month,
                            int day, int hour);

void geoCentricToDipole(MagSph& magCoords, GeoSph geoCoords,
                        Eigen::Matrix3d magneticRotationMatrix);

// Coordinate Conversions
GeoGeo convert(const GeoSph& sph);
MagGeo convert(const MagSph& sph);
GeoSph convert(const GeoGeo& geo);
MagSph convert(const MagGeo& geo);

// Solar Calculations
void computeSunPosition(double& alpha, double& declination, double ttTime);

double computeSolarZenith(double utTime, double ttTime, GeoGeo coords);

void computeSubSolar(GeoGeo& subSolar, double utTime, double ttTime);

double computeMagneticLocalTime(MagGeo subsolar, MagGeo observer);

double fourierSeries(nlohmann::json coefficients, double mlt);

double epsteinFunction(double h, double r, double h0, double S1, double S2);

void computeDipoleRotationMatrix(Eigen::Matrix3d& globalToMagRotation);
