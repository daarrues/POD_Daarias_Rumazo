//------------------------------------------------------------------------------
//                                    R_z
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file R_z.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/23
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero R_z.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "MatLabUtils/MatLabUtils.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void R_z(double angle, double rotmat[3][3])
//------------------------------------------------------------------------------
/**
 * Return a rotation matrix (of z).
 *
 * @param <angle> angle of rotation (rad)
 *
 * @return <rotmat> rotation matrix (3x3)
 */
//------------------------------------------------------------------------------
void R_z(double angle, double rotmat[3][3])
{
  double C = cos(angle);
  double S = sin(angle);

  zeros(rotmat[0],3);
  zeros(rotmat[1],3);
  zeros(rotmat[2],3);

  rotmat[0][0] =      C;  rotmat[0][1] =   S;  rotmat[0][2] = 0.0;
  rotmat[1][0] = -1.0*S;  rotmat[1][1] =   C;  rotmat[1][2] = 0.0;
  rotmat[2][0] =    0.0;  rotmat[2][1] = 0.0;  rotmat[2][2] = 1.0;
}
//------------------------------------------------------------------------------
