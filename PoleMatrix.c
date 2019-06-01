//------------------------------------------------------------------------------
//                                 PoleMatrix
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file PoleMatrix.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero PoleMatrix.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "PoleMatrix.h"
#include "R_x.h"
#include "R_y.h"
#include "MatLabUtils/MatLabUtils.h"

//------------------------------------------------------------------------------
//  void PoleMatrix(double xp, double yp, double PoleMat[3][3])
//------------------------------------------------------------------------------
/**
 * Transformation from pseudo Earth-fixed to Earth-fixed coordinates for
 * a given date.
 *
 * @param <xp> Pole coordinate (in).
 * @param <yp> Pole coordinate (in).
 * @param <PoleMat> Pole matrix (out).
 */
//------------------------------------------------------------------------------
void PoleMatrix(double xp, double yp, double PoleMat[3][3])
{
  double R_yResult[3][3];
  double R_xResult[3][3];

  R_y(-xp, R_yResult);
  R_x(-yp, R_xResult);
  prodMatr(R_yResult, R_xResult, PoleMat);
}
//------------------------------------------------------------------------------
