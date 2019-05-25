//------------------------------------------------------------------------------
//                                   gast
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file gast.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero gast.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "EOPDATA.h"

#include "SAT_Const.h"
#include "IERS.h"
#include "timediff.h"
#include "gmst.h"
#include "EqnEquinox.h"
#include <math.h>

//------------------------------------------------------------------------------
//  double gast(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * Greenwich Apparent Sidereal Time.
 *
 * @param <Mjd_UT1> Modified Julian Date UT1.
 *
 * @return GAST in [rad].
 */
//------------------------------------------------------------------------------
double gast(double Mjd_UT1)
{
  double *eopdata[13];
  leerFichero(eopdata);

  double UT1_UTC;
  double TAI_UTC;
  double x_pole;
  double y_pole;
  double ddpsi;
  double ddeps;

  IERS(eopdata, 20026, Mjd_UT1, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole,
       &ddpsi, &ddeps);

  double UT1_TAI;
  double UTC_GPS;
  double UT1_GPS;
  double TT_UTC;
  double GPS_UTC;

  timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

  double Mjd_UTC = Mjd_UT1 - UT1_UTC/86400;
  double Mjd_TT = Mjd_UTC + TT_UTC/86400;

  return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_TT), pi2);
}
