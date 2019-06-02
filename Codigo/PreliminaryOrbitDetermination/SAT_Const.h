//------------------------------------------------------------------------------
//                                SAT_Const
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file SAT_Const.h
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/09
 *
 * Este fichero contiene la declaración de constantes matemáticas
 * y astronómicas.
 * Reimplementación de SAT_Const.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#ifndef SATCONST_H
#define SATCONST_H

// Mathematical constants
const double pi2;
const double Rad;
const double Deg;
const double Arcs;

// General
const double MJD_J2000;
const double T_B1950;
const double c_light;
const double AU;

// Physical parameters of the Earth, Sun and Moon

// Equatorial radius and flattening
const double R_Earth;
const double f_Earth;
const double R_Sun;
const double R_Moon;

// Earth rotation (derivative of GMST at J2000;
// differs from inertial period by precession)
const double omega_Earth;

// Gravitational coefficients
const double GM_Earth;
const double GM_Sun;
const double GM_Moon;
const double GM_Mercury;
const double GM_Venus;
const double GM_Mars;
const double GM_Jupiter;
const double GM_Saturn;
const double GM_Uranus;
const double GM_Neptune;
const double GM_Pluto;

// Solar radiation pressure at 1 AU
const double P_Sol;

#endif // SATCONST_H
