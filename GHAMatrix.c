//------------------------------------------------------------------------------
//                                 GHAMatrix
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file GHAMatrix.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero GHAMatrix.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "R_z.h"
#include "gast.h"

//------------------------------------------------------------------------------
//  void GHAMatrix(double Mjd_UT1, double GHAmat[3][3])
//------------------------------------------------------------------------------
/**
 * Transformation from true equator and equinox to Earth equator and Greenwich
 * meridian system.
 *
 * @param <Mjd_UT1> Modified Julian Date UT1 (in).
 * @param <GHAmat> Greenwich Hour Angle matrix (out).
 */
//------------------------------------------------------------------------------
void GHAMatrix(double Mjd_UT1, double GHAmat[3][3])
{
  R_z(gast(Mjd_UT1), GHAmat);
}
