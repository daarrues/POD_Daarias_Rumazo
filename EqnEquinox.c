//------------------------------------------------------------------------------
//                                 EqnEquinox
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file EqnEquinox.h
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/16
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero EqnEquinox.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "NutAngles.h"
#include "MeanObliquity.h"
#include <math.h>

//------------------------------------------------------------------------------
//  double EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computation of the equation of the equinoxes.
 *
 * @param <Mjd_TT> Modified Julian Date (Terrestrial Time).
 *
 * @return Equation of the equinoxes.
 */
//------------------------------------------------------------------------------
double EqnEquinox(double Mjd_TT)
{
  double dpsi;
  double deps;

  // Nutation in longitude and obliquity
  NutAngles(Mjd_TT, &dpsi, &deps);

  // Equation of the equinoxes
  return dpsi * cos(MeanObliquity(Mjd_TT));
}
