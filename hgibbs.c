//------------------------------------------------------------------------------
//                              Hgibbs
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file hgibbs.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/01
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero hgibbs.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "angl.h"
#include "hgibbs.h"
#include "MatLabUtils/MatLabUtils.h"
#include "SAT_Const.h"
#include <math.h>
#include <string.h>

//------------------------------------------------------------------------------
//  void gibbs(double r1[], double r2[], double r3[], double MJD1, double MJD2,
//             double MJD3, double v2[], double *theta, double *theta1,
//             double *copa, char error[])
//------------------------------------------------------------------------------
/*
 *  hgibbs: implements the herrick-gibbs approximation for orbit
 *          determination, and finds the middle velocity vector for the 3
 *          given position vectors.
 *
 * Inputs:
 * @param <r1> - ijk position vector #1 (m)
 * @param <r2> - ijk position vector #2 (m)
 * @param <r3> - ijk position vector #3 (m)
 * @param <MJD1> - julian date of 1st sighting (days since 4713 bc)
 * @param <MJD2> - julian date of 2nd sighting (days since 4713 bc)
 * @param <MJD3> - julian date of 3rd sighting (days since 4713 bc)
 *
 * Outputs:
 * @return <v2> - ijk velocity vector for r2 (m/s)
 * @return <theta>- angl between vectors (rad)
 * @return <error> - flag indicating success
 */
//------------------------------------------------------------------------------
void gibbs(double r1[], double r2[], double r3[], double MJD1, double MJD2,
           double MJD3, double v2[], double *theta, double *theta1,
           double *copa, char error[])
{
  strcpy(error, "          ok");
  *theta = 0.0;
  *theta1 = 0.0;
  double magr1 = norm(r1);
  double magr2 = norm(r2);
  double magr3 = norm(r3);
  zeros(v2, 3);

  double tolangle = 0.01745329251994;
  double dt21 = (MJD2-MJD1)*86400;
  double dt31= (MJD3-MJD1)*86400;
  double dt32= (MJD3-MJD2)*86400;

  double p[3], pn[3], r1n[3];
  cross(r2, r3, p);
  unit(p, pn);
  unit(r1, r1n);
  *copa = asin(dot(pn, r1n));

  if ( fabs(dot(r1n, pn)) > 0.017452406 )
  {
    strcpy(error, "not coplanar");
  }

  *theta  = angl(r1, r2);
  *theta1 = angl(r2, r3);

  if ((theta > tolangle) || (theta1 > tolangle))
  {
    strcpy(error, "   angl > 1ø");
  }

  double term1 = -dt32*(1/(dt21*dt31) + GM_Earth/(12*magr1*magr1*magr1));
  double term2 = (dt32-dt21)*(1/(dt21*dt32) + GM_Earth/(12*magr2*magr2*magr2));
  double term3 =  dt21*(1/(dt32*dt31) + GM_Earth/(12*magr3*magr3*magr3));

  for(int i = 0; i < 3; i++)
  {
    v2[i] = term1*r1[i] + term2*r2[i] + term3*r3[i];
  }
}
//------------------------------------------------------------------------------
