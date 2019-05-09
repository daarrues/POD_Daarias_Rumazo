//------------------------------------------------------------------------------
//                              lambert_gooding
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file lambert_gooding.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/09
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero lambert_gooding.m (M. Mahooti)
 */
//------------------------------------------------------------------------------

#include "MatLabUtils/MatLabUtils.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void lambert_gooding(double v[3])
//------------------------------------------------------------------------------
/**
* Lambert's problem using Gooding's method
*
* @param <r1>            first cartesian position [km]
* @param <r2>            second cartesian position [km]
* @param <tof>           time of flight [sec]
* @param <mu>            gravity parameter [km^3/s^2]
* @param <long_way>      when true, do "long way" (>pi) transfers
* @param <multi_revs>    maximum number of multi-rev solutions to compute
* @param <v1>           return, vector containing 3d arrays with the cartesian components
*                of the velocities at r1
* @param <v2>            return, vector containing 3d arrays with the cartesian components
*                of the velocities at r2
*
* @note References:
*  1. R. H, Gooding. "[A procedure for the solution of Lambert's orbital
*     boundary-value problem](http://adsabs.harvard.edu/abs/1990CeMDA..48..145G)"
*     Celestial Mechanics and Dynamical Astronomy,
*     vol. 48, no. 2, 1990, p. 145-165.
*  2. A. Klumpp, "Performance Comparision of Lambert and Kepler Algorithms",
*     JPL Interoffice Memorandum, 314.1-0426-ARK, Jan 2, 1991.
*     [Zip](http://derastrodynamics.com/docs/lambert_papers_v1.zip)
*
*/
//------------------------------------------------------------------------------
void lambert_gooding(double r1, double r2, double tof, int long_way,
                     int multi_revs, double v1[][], double v2[][])
{

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//  void lambert_gooding(double v[3])
//------------------------------------------------------------------------------
/**
* 8th root function, used by xlamb
*
* @param <x>            number to find the 8th root to
*
* @return the 8th root of x
*/
//------------------------------------------------------------------------------
double d8rt (double x)
{
  // 8th root function, used by xlamb
  return sqrt(sqrt(sqrt(x)));
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  void lambert_gooding(double v[3])
//------------------------------------------------------------------------------
/**
* 8th root function, used by xlamb
*
* @param <x>            number to find the 8th root to
*
* @return the 8th root of x
*/
//------------------------------------------------------------------------------
void xlamb(m,q,qsqfm1,tin)
//------------------------------------------------------------------------------
