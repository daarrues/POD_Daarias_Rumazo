//------------------------------------------------------------------------------
//                                 PrecMatrix
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file PrecMatrix.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero PrecMatrix.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "R_y.h"
#include "R_z.h"
#include "SAT_Const.h"
#include "MatLabUtils/MatLabUtils.h"

//------------------------------------------------------------------------------
//  void PrecMatrix(double Mjd_1, double Mjd_2, double PrecMat[3][3])
//------------------------------------------------------------------------------
/**
 * Precession transformation of equatorial coordinates.
 *
 * @param <Mjd_1> Epoch given (Modified Julian Date TT) (in).
 * @param <Mjd_2> Epoch to precess to (Modified Julian Date TT) (in).
 * @param <PrecMat> Precession transformation matrix (out).
 */
//------------------------------------------------------------------------------
void PrecMatrix(double Mjd_1, double Mjd_2, double PrecMat[3][3])
{
  double T  = (Mjd_1 - MJD_J2000)/36525;
  double dT = (Mjd_2 - Mjd_1)/36525;

  // Precession angles
  double zeta  = ((2306.2181 + (1.39656 - 0.000139*T)*T) +
                 ((0.30188 - 0.000344*T) + 0.017998*dT)*dT)*dT/Arcs;
  double z     = zeta + ((0.79280 + 0.000411*T) + 0.000205*dT)*dT*dT/Arcs;
  double theta = ((2004.3109 - (0.85330 + 0.000217*T)*T) -
                 ((0.42665 + 0.000217*T) + 0.041833*dT)*dT)*dT/Arcs;

  // Precession matrix
  double R_zResult1[3][3];
  double R_yResult[3][3];
  double R_zResult2[3][3];
  double prodResult[3][3];

  R_z(-z, R_zResult1);
  R_y(theta, R_yResult);
  R_z(-zeta, R_zResult2);
  prodMatr(R_zResult1, R_yResult, prodResult);
  prodMatr(prodResult, R_zResult2, PrecMat);
}
//------------------------------------------------------------------------------
