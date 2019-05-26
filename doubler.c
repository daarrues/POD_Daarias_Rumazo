//------------------------------------------------------------------------------
//                                    doubler
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file doubler.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/26
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero doubler.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "MatLabUtils/MatLabUtils.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void doubler(double cc1, double cc2, double magrsite1, double magrsite2,
//           double magr1in, double magr2in,
//           double los1[3], double los2[3], double los3[3],
//           double rsite1[3], double rsite2[3], double rsite3[3],
//           double t1, double t3, char direct,
//           double r2[], double r3[], double *f1, double *f2, double *q1,
//           double *magr1, double *magr2, double *a, double *deltae32)
//------------------------------------------------------------------------------
/**
 * this rountine accomplishes the iteration work for the double-r angles
 * only routine
 *
 * Inputs:
 * @param <cc1> cc1
 * @param <cc2> cc2
 * @param <magrsite1> magrsite1
 * @param <magrsite2> magrsite2
 * @param <magr1in> magr1in
 * @param <magr2in> magr2in
 * @param <los1> los1
 * @param <los2> los2
 * @param <los3> los3
 * @param <rsite1> rsite1
 * @param <rsite2> rsite2
 * @param <rsite3> rsite3
 * @param <t1> t1
 * @param <t3> t3
 * @param <direct> direct
 *
 * Output:
 * @param <r2> r2
 * @param <r3> r3
 * @param <f1> f1
 * @param <f2> f2
 * @param <q1> q1
 * @param <magr1> magr1
 * @param <magr2> magr2
 * @param <a> a
 * @param <deltae32> deltae32
 */
//------------------------------------------------------------------------------
void doubler(double cc1, double cc2, double magrsite1, double magrsite2,
             double magr1in, double magr2in,
             double los1[3], double los2[3], double los3[3],
             double rsite1[3], double rsite2[3], double rsite3[3],
             double t1, double t3, char direct,
             double r2[3], double r3[3], double *f1, double *f2, double *q1,
             double *magr1, double *magr2, double *a, double *deltae32)
{
  double mu = 398600.4418e+9;

  double rho1 = (-cc1+sqrt(pow(cc1,2)-4*(pow(magrsite1,2)-pow(magr1in,2))))/2;
  double rho2 = (-cc2+sqrt(pow(cc2,2)-4*(pow(magrsite2,2)-pow(magr2in,2))))/2;

  double r1[3];
  for (int i = 0; i < 3; i++) {
    r1[i] = rho1 * los1[i] + rsite1[i];
    r2[i] = rho2 * los2[i] + rsite2[i];
  }
  *magr1 = norm(r1);
  *magr2 = norm(r2);

  double w[3];
  double crossR1R2[3];
  cross(r1, r2, crossR1R2);
  for(int i = 0; i < 3; i++)
  {
    if(direct == 'y')
    {
      w[i] = crossR1R2[i] / (*magr1 * *magr2);
    }
    else
    {
      w[i] = -crossR1R2[i] / (*magr1 * *magr2);
    }
  }

  // change to negative sign
  double rho3 = -dot(rsite3, w) / dot(los3, w);
  for (int i = 0; i < 3; i++) {
    r3[i] = rho3 * los3[i] + rsite3[i];
  }
  double magr3 = norm(r3);

  double cosdv21 = dot(r2, r1) / (*magr2 * *magr1);
  double sindv21 = sqrt(1 - pow(cosdv21,2));
  double dv21 = atan2(sindv21, cosdv21);

  double cosdv31 = dot(r3, r1) / (magr3 * *magr1);
  double sindv31 = sqrt(1 - pow(cosdv31,2));
  double dv31 = atan2(sindv31, cosdv31);

  double cosdv32 = dot(r3, r2) / (magr3 * *magr2);
  double sindv32 = sqrt(1 - pow(cosdv32,2));

  double c1;
  double c3;
  double p;
  if(dv31 > M_PI)
  {
    c1 = (*magr2 * sindv32) / (*magr1 * sindv31);
    c3 = (*magr2 * sindv21) / (magr3 * sindv31);
    p = (c1 * *magr1 + c3 * magr3 - *magr2) / (c1 + c3 - 1);
  }
  else
  {
    c1 = (*magr1 * sindv31) / (*magr2 * sindv32);
    c3 = (*magr1 * sindv21) / (magr3 * sindv32);
    p = (c3 * magr3 - c1 * *magr2 + *magr1) / (-c1 + c3 + 1);
  }

  double ecosv1 = p / *magr1 - 1;
  double ecosv2 = p / *magr2 - 1;
  double ecosv3 = p / magr3 - 1;

  double esinv2;
  if(fabs(dv21 - M_PI) < 10e-12)
  {
    esinv2 = (-cosdv21*ecosv2 + ecosv1) / sindv21;
  }
  else
  {
    esinv2 = (cosdv32*ecosv2 - ecosv3) / sindv31;
  }

  double e = sqrt(pow(ecosv2,2) + pow(esinv2,2));
  *a = p / (1 - pow(e,2));

  double n;
  double s;
  double c;
  double deltam32;
  double deltam12;
  if(pow(e,2) < 1)
  {
    n = sqrt(mu / pow(*a,3));
    s = *magr2 / p * sqrt(1 - pow(e,2)) * esinv2;
    c = *magr2 / p * (pow(e,2) + ecosv2);

    double sinde32 = magr3/sqrt(*a*p)*sindv32 - magr3/p*(1 - cosdv32)*s;
    double cosde32 = 1 - *magr2*magr3/(*a*p)*(1 - cosdv32);
    *deltae32 = atan2(sinde32, cosde32);

    double sinde21 = *magr1/sqrt(*a*p)*sindv21 + *magr1/p*(1 - cosdv21)*s;
    double cosde21 = 1 - *magr2**magr1/(*a*p)*(1 - cosdv21);
    double deltae21 = atan2(sinde21, cosde21);

    deltam32 = *deltae32 + 2*s*pow(sin(*deltae32/2),2) - c*sin(*deltae32);
    deltam12 = -deltae21 + 2*s*pow(sin(deltae21/2),2) + c*sin(deltae21);
  }
  else
  {
    n = sqrt(mu / (-pow(*a,3)));
    s = *magr2 / p * sqrt(pow(e,2) - 1) * esinv2;
    c = *magr2 / p * (pow(e,2) + ecosv2);

    double sindh32 = magr3/sqrt(-*a*p)*sindv32 - magr3/p*(1 - cosdv32)*s;
    double sindh21 = *magr1/sqrt(-*a*p)*sindv21 + *magr1/p*(1 - cosdv21)*s;

    double deltah32 = log(sindh32 + sqrt(pow(sindh32,2) + 1));
    double deltah21 = log(sindh21 + sqrt(pow(sindh21,2) + 1));

    deltam32 = -deltah32 + 2*s*pow(sinh(deltah32/2),2) + c*sinh(deltah32);
    deltam12 = deltah21 + 2*s*pow(sinh(deltah21/2),2) - c*sinh(deltah21);
    *deltae32 = deltah32;
  }

  *f1 = t1 - deltam12/n;
  *f2 = t3 - deltam32/n;

  *q1 = sqrt(pow(*f1,2) + pow(*f2,2));
}
//------------------------------------------------------------------------------
