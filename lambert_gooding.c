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
//  void lambert_gooding(double r1[], double r2[], double tof, double long_way,
//                       double multi_revs, double v1[], double v2[])
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
void lambert_gooding(double r1[], double r2[], double tof, double long_way,
                     double multi_revs, double v1[], double v2[])
{
  // temp arrays to hold all the solutions:
  // they will be packed into the output arrays
  // logical,dimension(2*multi_revs+1) :: solution_exists
  // real(wp),dimension(3,1+2*multi_revs) :: all_vt1, all_vt2

  double r1mag = norm(r1);
  double r2mag = norm(r2);

  if ( r1mag==0.0 || r2mag==0.0 || mu<=0.0 || tof<=0.0 )
  {
      print("Error in solve_lambert_gooding: invalid input\n");
      return;
  }

  // initialize:
  int solution_exists = 0;
  double dr       = r1mag - r2mag;
  double r1r2     = r1mag*r2mag;
  double r1hat[3] = unit(r1);
  double r2hat[3] = unit(r2);
  double r1xr2    = cross(r1,r2);
  if (all(r1xr2==0.0)) // the vectors are parallel,
  {                    // so the transfer plane is undefined
      r1xr2 = {0.0, 0.0, 1.0}; // degenerate conic...choose the x-y plane
  }
  double r1xr2_hat[3] = unit(r1xr2);

  // a trick to make sure argument is between [-1 and 1]:
  double pa = acos(max(-1.0,min(1.0,dot(r1hat,r2hat))));

  int num_revs;
  double ta, rho[3], etai[3], etaf[3]
  for (int i = 0; i < multi_revs; i++)
  {
      num_revs = i; // number of complete revs for this case

      // transfer angle and normal vector:
      if (long_way) // greater than pi
      {
          ta    =  num_revs * 2*M_PI + (2*M_PI - pa);
          rho   = -r1xr2_hat;
      }
      else // less than pi
      {
          ta    = num_revs * 2*M_PI + pa;
          rho   = r1xr2_hat;
      }
      etai = cross(rho,r1hat);
      etaf = cross(rho,r2hat);

      // Gooding routine:
      double n, vri[2], vti[2], vrf[2], vtf[2];
      vlamb(mu, r1mag, r2mag, ta, tof, n, vri, vti, vrf, vtf);
      double vt1[3][n], vt2[3][n];

      switch (n) // number of solutions
      {
        case 1:
          for(int j = 0; j < 3; j++)
          {
            vt1[j][0] = vri[0]*r1hat[j] + vti[0]*etai[j];
            vt2[j][0] = vrf[0]*r2hat[j] + vtf[0]*etaf[j];
          }
          break;
        case 2:
          for(int j = 0; j < 3; j++)
          {
            vt1[j][0] = vri[0]*r1hat[j] + vti[0]*etai[j];
            vt2[j][0] = vrf[0]*r2hat[j] + vtf[0]*etaf[j];
            vt1[j][1] = vri[1]*r1hat[j] + vti[1]*etai[j];
            vt2[j][1] = vrf[1]*r2hat[j] + vtf[1]*etaf[j];
          }
          break;
      }

      if (i == 0 && n == 1) // there can be only one solution
      {
          all_vt1(:,1) = vt1(:,1);
          all_vt2(:,1) = vt2(:,1);
          solution_exists(1) = true;
      }
      else
      {
          switch (n)
          {
              case 1:
                  all_vt1(:,2*i)         = vt1(:,1);
                  all_vt2(:,2*i)         = vt2(:,1);
                  solution_exists(2*i)   = true;
                  break;
              case 2:
                  all_vt1(:,2*i)         = vt1(:,1);
                  all_vt2(:,2*i)         = vt2(:,1);
                  solution_exists(2*i)   = true;
                  all_vt1(:,2*i+1)       = vt1(:,2);
                  all_vt2(:,2*i+1)       = vt2(:,2);
                  solution_exists(2*i+1) = true;
                  break;
          }
      }
  }

  // return all the solutions:
  n_solutions = length(solution_exists);

  v1 = zeros(3,n_solutions);
  v2 = zeros(3,n_solutions);

  k=0;
  for i=1:size(solution_exists)
  {
      if (solution_exists(i))
      {
          k=k+1;
          v1(:,k) = all_vt1(:,i);
          v2(:,k) = all_vt2(:,i);
      }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  void vlamb(double gm, double r1, double r2, double th, double tdelt,
//             double n, double vri[], double vti[], double vrf[], double vtf[])
//------------------------------------------------------------------------------
/**
 * Gooding support routine
 * @note this contains the modification from [2]
 *
 * @param <x>            number to find the 8th root to
 *
 * @return the 8th root of x
 */
//------------------------------------------------------------------------------
void vlamb(double gm, double r1, double r2, double th, double tdelt, double n,
           double vri[], double vti[], double vrf[], double vtf[])
{

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  void tlamb(double m, double q, double qsqfm1, double x, double n, double t,
//             double dt, double d2t, double d3t)
//------------------------------------------------------------------------------
/**
 * Gooding support routine
 *
 * @param <x>            number to find the 8th root to
 *
 * @return the 8th root of x
 */
//------------------------------------------------------------------------------
void tlamb(double m, double q, double qsqfm1, double x, double n, double t,
           double dt, double d2t, double d3t)
{
  double y, z, qx, a, b, aa, bb, g, f, fg1, term, fg1sq, twoi1, told, qz, qz2,
         u0i, u1i, u2i, u3i, tq, i, tqsum, ttmold, p, tterm, tqterm;

  // Gooding support routine
  double sw = 0.4;
  t = 0;

  double lm1 = (n == -1);
  double l1 = (n >= 1);
  double l2 = (n >= 2);
  double l3 = (n == 3);
  double qsq = q*q;
  double xsq = x*x;
  double u = (1.0 - x)*(1.0 + x);

  if (!lm1)
  {
      // (needed if series, and otherwise useful when z = 0)
      dt = 0.0;
      d2t = 0.0;
      d3t = 0.0;
  }

  if (lm1 || m > 0 || x < 0.0 || fabs(u) > sw)
  {
    // direct computation (not series)
    y = sqrt(fabs(u));
    z = sqrt(qsqfm1 + qsq*xsq);
    qx = q*x;

    if (qx <= 0.0)
    {
        a = z - qx;
        b = q*z - x;
    }
    if (qx < 0.0 && lm1)
    {
        aa = qsqfm1/a;
        bb = qsqfm1*(qsq*u - xsq)/b;
    }
    if (qx == 0.0 && lm1 || qx > 0.0)
    {
        aa = z + qx;
        bb = q*z + x;
    }
    if (qx > 0.0)
    {
        a = qsqfm1/aa;
        b = qsqfm1*(qsq*u - xsq)/bb;
    }
    if (!lm1)
    {
      if (qx*u >= 0.0)
      {
        g = x*z + q*u;
      }
      else
      {
          g = (xsq - qsq*u)/(x*z - q*u);
      }
      f = a*y;
      if (x <= 1.0)
      {
        t = m*M_PI + atan2(f, g);
      }
      else
      {
        if (f > sw)
        {
            t = log(f + g);
        }
        else
        {
          fg1 = f/(g + 1.0);
          term = 2.0*fg1;
          fg1sq = fg1*fg1;
          t = term;
          twoi1 = 1.0;
          while(1)
          {
            twoi1 = twoi1 + 2.0;
            term = term*fg1sq;
            told = t;
            t = t + term/twoi1;
            if (t != told)
            {
              // cycle
              break;
            }
          }// (continue looping for inverse tanh)
        }
      }
      t = 2.0*(t/y + b)/u;
      if (l1 && z != 0.0)
      {
        qz = q/z;
        qz2 = qz*qz;
        qz = qz*qz2;
        dt = (3.0*x*t - 4.0*(a + qx*qsqfm1)/z)/u;
        if (l2)
        {
            d2t = (3.0*t + 5.0*x*dt + 4.0*qz*qsqfm1)/u;
        }
        if (l3)
        {
            d3t = (8.0*dt + 7.0*x*d2t - 12.0*qz*qz2*x*qsqfm1)/u;
        }
      }
    }
    else
    {
      dt = b;
      d2t = bb;
      d3t = aa;
    }
  }
  else
  {
      // compute by series
      u0i = 1.0;
      if (l1)
      {
        u1i = 1.0;
      }
      if (l2)
      {
        u2i = 1.0;
      }
      if (l3)
      {
        u3i = 1.0;
      }
      term = 4.0;
      tq = q*qsqfm1;
      i = 0;
      if (q<0.5)
      {
        tqsum = 1.0 - q*qsq;
      }
      if (q>=0.5)
      {
        tqsum = (1.0/(1.0 + q) + q)*qsqfm1;
      }
      ttmold = term/3.0;
      t = ttmold*tqsum;
      while(1)
      {
        i = i + 1;
        p = i;
        u0i = u0i*u;
        if (l1 && i>1)
        {
            u1i = u1i*u;
        }
        if (l2 && i>2)
        {
            u2i = u2i*u;
        }
        if (l3 && i>3)
        {
            u3i = u3i*u;
        }
        term = term*(p - 0.5)/p;
        tq = tq*qsq;
        tqsum = tqsum + tq;
        told = t;
        tterm = term/(2.0*p + 3.0);
        tqterm = tterm*tqsum;
        t = t - u0i*((1.5*p + 0.25)*tqterm/(p*p - 0.25)-ttmold*tq);
        ttmold = tterm;
        tqterm = tqterm*p;
        if (l1)
        {
          dt = dt + tqterm*u1i;
        }
        if (l2)
        {
          d2t = d2t + tqterm*u2i*(p - 1.0);
        }
        if (l3)
        {
          d3t = d3t + tqterm*u3i*(p - 1.0)*(p - 2.0);
        }
        if (i < n || t != told)
        {
          // cycle
          break;
        }
      }
    if (l3)
    {
      d3t = 8.0*x*(1.5*d2t - xsq*d3t);
    }
    if (l2)
    {
      d2t = 2.0*(2.0*xsq*d2t - dt);
    }
    if (l1)
    {
      dt = -2.0*x*dt;
    }
    t = t/xsq;
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  double d8rt (double x)
//------------------------------------------------------------------------------
/**
 * 8th root function, used by xlamb
 *
 * @param <x> number to find the 8th root to
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
//  void xlamb(double v[3])
//------------------------------------------------------------------------------
/**
 * 8th root function, used by xlamb
 *
 * @param <x>            number to find the 8th root to
 *
 * @return the 8th root of x
 */
//------------------------------------------------------------------------------
void xlamb(double m, double q, double qsqfm1, double tin, double n, double x,
           double xpl)
{

}
//------------------------------------------------------------------------------