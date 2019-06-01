//------------------------------------------------------------------------------
//                                  Frac
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file Frac.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/22
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero Frac.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "Frac.h"
#include <math.h>

//------------------------------------------------------------------------------
//  double Frac(double n)
//------------------------------------------------------------------------------
/**
 * Fractional part of a number (y=x-[x])
 *
 * @param <n> A real number.
 *
 * @return Fractional part of n.
 */
//------------------------------------------------------------------------------
double Frac(double n)
{
  return n - floor(n);
}
//------------------------------------------------------------------------------
