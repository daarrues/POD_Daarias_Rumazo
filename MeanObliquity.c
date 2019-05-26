//------------------------------------------------------------------------------
//                              MeanObliquity
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file MeanObliquity.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/15
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero MeanObliquity.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "SAT_Const.h"

//------------------------------------------------------------------------------
//  double MeanObliquity(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computes the mean obliquity of the ecliptic.
 *
 * @param <Mjd_TT> Modified Julian Date (Terrestrial Time).
 *
 * @return Mean obliquity of the ecliptic.
 */
//------------------------------------------------------------------------------
double MeanObliquity(double Mjd_TT)
{
  double T = (Mjd_TT - MJD_J2000) / 36525;
  return Rad*(23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600);
}
//------------------------------------------------------------------------------
