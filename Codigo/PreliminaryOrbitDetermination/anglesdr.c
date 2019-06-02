//------------------------------------------------------------------------------
//                                    anglesdr
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file anglesdr.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/01
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero anglesdr.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "../MatLabUtils/MatLabUtils.h"
#include "anglesdr.h"
#include "doubler.h"
#include "lambert_gooding.h"
#include "SAT_Const.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1,
//         double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3,
//         double rsite1[3], double rsite2[3], double rsite3[3],
//         double r2[3], double v2[3])
//------------------------------------------------------------------------------
/**
 * solves the problem of orbit determination using three optical sightings.
 *
 * Inputs:
 * @param <rtasc1> right ascension at t1 (rad)
 * @param <rtasc2> right ascension at t2 (rad)
 * @param <rtasc3> right ascension at t3 (rad)
 * @param <decl1> declination at t1 (rad)
 * @param <decl2> declination at t2 (rad)
 * @param <decl3> declination at t3 (rad)
 * @param <Mjd1> Modified julian date of t1
 * @param <Mjd2> Modified julian date of t2
 * @param <Mjd3> Modified julian date of t3
 * @param <rsite1> ijk site1 position vector (m)
 * @param <rsite2> ijk site2 position vector (m)
 * @param <rsite3> ijk site3 position vector (m)
 *
 * Output:
 * @param <r> ijk position vector at t2 (m)
 * @param <v> ijk velocity vector at t2 (m/s)
 */
//------------------------------------------------------------------------------
void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1,
           double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3,
           double rsite1[3], double rsite2[3], double rsite3[3],
           double r2[3], double v2[3])
{
  double magr1in = 2.01 * R_Earth;
  double magr2in = 2.11 * R_Earth;
  char direct = 'y';

  double tol    = 1e-8 * R_Earth;
  double pctchg = 5e-6;

  double t1 = (Mjd1 - Mjd2)*86400;
  double t3 = (Mjd3 - Mjd2)*86400;

  double los1[3] = {cos(decl1)*cos(rtasc1), cos(decl1)*sin(rtasc1), sin(decl1)};
  double los2[3] = {cos(decl2)*cos(rtasc2), cos(decl2)*sin(rtasc2), sin(decl2)};
  double los3[3] = {cos(decl3)*cos(rtasc3), cos(decl3)*sin(rtasc3), sin(decl3)};

  double magr1old  = 99999e3;
  double magr2old  = 99999e3;
  double magrsite1 = norm(rsite1);
  double magrsite2 = norm(rsite2);
  double magrsite3 = norm(rsite3);

  double cc1 = 2*dot(los1, rsite1);
  double cc2 = 2*dot(los2, rsite2);


  double r3[3];
  double f1;
  double f2;
  double q1;
  double magr1;
  double magr2;
  double a;
  double deltae32;
  double f;
  double g;
  int ll = 0;
  while(fabs(magr1in-magr1old) > tol && fabs(magr2in-magr2old) > tol && ll <= 3)
  {
    ll++;
    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
      los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct,
      r2, r3, &f1, &f2, &q1, &magr1, &magr2, &a, &deltae32);

    f = 1 - a/magr2 * (1 - cos(deltae32));
    g = t3 - sqrt(pow(a,3)/GM_Earth) * (deltae32 - sin(deltae32));
    for (int i = 0; i < 3; i++)
    {
      v2[i] = (r3[i] - f*r2[i]) / g;
    }

    double magr1o = magr1in;
    magr1in = (1 + pctchg) * magr1in;
    double deltar1 = pctchg * magr1in;
    double f1delr1;
    double f2delr1;
    double q2;
    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
      los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct,
      r2, r3, &f1delr1, &f2delr1, &q2, &magr1, &magr2, &a, &deltae32);

    double pf1pr1 = (f1delr1 - f1) / deltar1;
    double pf2pr1 = (f2delr1 - f2) / deltar1;

    magr1in = magr1o;
    deltar1 = pctchg * magr1in;
    double magr2o = magr2in;
    magr2in = (1 + pctchg) * magr2in;
    double deltar2 = pctchg * magr2in;
    double f1delr2;
    double f2delr2;
    double q3;
    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
      los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct,
      r2, r3, &f1delr2, &f2delr2, &q3, &magr1, &magr2, &a, &deltae32);

    double pf1pr2 = (f1delr2 - f1) / deltar2;
    double pf2pr2 = (f2delr2 - f2) / deltar2;

    magr2in = magr2o;
    deltar2 = pctchg * magr2in;

    double delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
    double delta1 = pf2pr2*f1 - pf1pr2*f2;
    double delta2 = pf1pr1*f2 - pf2pr1*f1;

    deltar1 = -delta1/delta;
    deltar2 = -delta2/delta;

    magr1old = magr1in;
    magr2old = magr2in;

    magr1in = magr1in + deltar1;
    magr2in = magr2in + deltar2;
  }

  doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
    los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct,
    r2, r3, &f1, &f2, &q1, &magr1, &magr2, &a, &deltae32);

  double v3[3];
  lambert_gooding(r2, r3, t3, GM_Earth, 0.0, 1.0, v2, v3);
  f = 1 - a/magr2 * (1 - cos(deltae32));
  g = t3 - sqrt(pow(a,3)/GM_Earth) * (deltae32 - sin(deltae32));
  for (int i = 0; i < 3; i++)
  {
    v2[i] = (r3[i] - f*r2[i]) / g;
  }
}
//------------------------------------------------------------------------------
