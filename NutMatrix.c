//------------------------------------------------------------------------------
//                                 NutMatrix
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file NutMatrix.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero NutMatrix.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"
#include "MatLabUtils/MatLabUtils.h"

//------------------------------------------------------------------------------
//  void NutMatrix(double Mjd_TT, double NutMat[3][3])
//------------------------------------------------------------------------------
/**
 * Transformation from mean to true equator and equinox.
 *
 * @param <Mjd_TT> Modified Julian Date (Terrestrial Time) (in).
 * @param <PoleMat> Nutation matrix (out).
 */
//------------------------------------------------------------------------------
void NutMatrix(double Mjd_TT, double NutMat[3][3])
{
  // Mean obliquity of the ecliptic
  double ep = MeanObliquity(Mjd_TT);

  // Nutation in longitude and obliquity
  double dpsi;
  double deps;
  NutAngles(Mjd_TT, &dpsi, &deps);

  // Transformation from mean to true equator and equinox
  double R_xResult1[3][3];
  double R_zResult[3][3];
  double R_xResult2[3][3];
  double prodResult[3][3];

  R_x(-ep-deps, R_xResult1);
  R_z(-dpsi, R_zResult);
  R_x(+ep, R_xResult2);
  prodMatr(R_xResult1, R_zResult, prodResult);
  prodMatr(prodResult, R_xResult2, NutMat);
}
//------------------------------------------------------------------------------
