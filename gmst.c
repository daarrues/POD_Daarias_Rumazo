//------------------------------------------------------------------------------
//                                  gmst
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file gmst.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/22
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero gmst.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "SAT_Const.h"
#include "Frac.h"
#include <math.h>

//------------------------------------------------------------------------------
//  double gmst(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * Greenwich Mean Sidereal Time.
 *
 * @param <Mjd_UT1> Modified Julian Date UT1.
 *
 * @return GMST in [rad].
 */
//------------------------------------------------------------------------------
double gmst(double Mjd_UT1)
{
  int Secs = 86400;                      // Seconds per day

  double Mjd_0 = floor(Mjd_UT1);
  double UT1   = Secs * (Mjd_UT1 - Mjd_0);      // [s]
  double T_0   = (Mjd_0 - MJD_J2000) / 36525;
  double T     = (Mjd_UT1 - MJD_J2000) / 36525;

  double gmst = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1
         + (0.093104 - 6.2e-6*T)*T*T;    // [s]

  return pi2*Frac(gmst/Secs);           // [rad], 0..2pi
}
//------------------------------------------------------------------------------
