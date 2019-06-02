//------------------------------------------------------------------------------
//                              Gibbs
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file gibbs.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/01
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero gibbs.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "../MatLabUtils/MatLabUtils.h"
#include "angl.h"
#include "gibbs.h"
#include "SAT_Const.h"
#include <math.h>
#include <string.h>

//------------------------------------------------------------------------------
//  void gibbs(double r1[], double r2[], double r3[], double v2[],
//             double *theta, double *theta1, double *copa, char error[])
//------------------------------------------------------------------------------
/*
 *  gibbs: performs the gibbs method of orbit determination. this method
 *         determines the velocity at the middle point of the 3 given
 *         position vectors.
 *
 * Inputs:
 * @param <r1> - ijk position vector #1 (m)
 * @param <r2> - ijk position vector #2 (m)
 * @param <r3> - ijk position vector #3 (m)
 *
 * Outputs:
 * @return <v2> - ijk velocity vector for r2 (m/s)
 * @return <theta>- angl between vectors (rad)
 * @return <error> - flag indicating success
 */
//------------------------------------------------------------------------------
void gibbs(double r1[], double r2[], double r3[], double v2[], double *theta,
           double *theta1, double *copa, char error[])
{
  double small = 0.00000001;
  *theta = 0.0;
  *theta1 = 0.0;
  zeros(v2, 3);
  strcpy(error, "          ok");

  double magr1 = norm(r1);
  double magr2 = norm(r2);
  double magr3 = norm(r3);

  double p[3], q[3], w[3], pn[3], r1n[3];
  cross(r2, r3, p);
  cross(r3, r1, q);
  cross(r1, r2, w);
  unit(p, pn);
  unit(r1, r1n);
  *copa = asin(dot(pn, r1n));

  if ( fabs(dot(r1n, pn)) > 0.017452406 )
  {
    strcpy(error, "not coplanar");
  }

  double d[3], n[3], dn[3], nn[3];
  for(int i = 0; i < 3; i++)
  {
    d[i] = p[i] + q[i] + w[i];
    n[i] = magr1*p[i] + magr2*q [i] + magr3*w[i];
  }
  double magd = norm(d);
  double magn = norm(n);
  unit(d, dn);
  unit(n, nn);

  // -------------------------------------------------------------
  // determine if  the orbit is possible. both d and n must be in
  // the same direction, and non-zero.
  // -------------------------------------------------------------
  if ((fabs(magd) < small) || (fabs(magn) < small) || (dot(nn, dn) < small))
  {
    strcpy(error, "impossible  ");
  }
  else
  {
    *theta  = angl(r1, r2);
    *theta1 = angl(r2, r3);

    // ----------- perform gibbs method to find v2 -----------
    double r1mr2 = magr1 - magr2;
    double r3mr1 = magr3 - magr1;
    double r2mr3 = magr2 - magr3;
    double s[3], b[3];
    cross(d, r2, b);
    double l  = sqrt(GM_Earth/(magd*magn));
    double tover2 = l/magr2;
    for(int i = 0; i < 3; i++)
    {
      s[i] = r1mr2*r3[i] + r3mr1*r2[i] + r2mr3*r1[i];
      v2[i] = tover2 * b[i] + l * s[i];
    }
  }
}
//------------------------------------------------------------------------------
