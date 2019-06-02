//------------------------------------------------------------------------------
//                                   rv2coe
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file rv2coe.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/01
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero rv2coe.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "../MatLabUtils/MatLabUtils.h"
#include "angl.h"
#include "newtonnu.h"
#include "rv2coe.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void rv2coe(double r[3], double v[3], double *p, double *a, double *ecc,
//           double *incl, double *omega, double *argp, double *nu, double *m,
//           double *arglat, double *truelon, double *lonper)
//------------------------------------------------------------------------------
/**
 * Finds the classical orbital elements given the geocentric equatorial
 * position and velocity vectors.
 *
 * Inputs:
 * @param <r> ijk position vector (m).
 * @param <v> ijk velocity vector (m/s).
 *
 * Outputs:
 * @param <p> semilatus rectum (m).
 * @param <a> semimajor axis (m).
 * @param <ecc> eccentricity.
 * @param <incl> inclination (0.0 to pi rad).
 * @param <omega> longitude of ascending node (0.0 to 2pi rad).
 * @param <argp> argument of perigee (0.0 to 2pi rad).
 * @param <nu> true anomaly (0.0 to 2pi rad).
 * @param <m> mean anomaly (0.0 to 2pi rad).
 * @param <arglat> argument of latitude (0.0 to 2pi rad).
 * @param <truelon> true longitude (0.0 to 2pi rad).
 * @param <lonper> longitude of periapsis (0.0 to 2pi rad).
 */
//------------------------------------------------------------------------------
void rv2coe(double r[3], double v[3], double *p, double *a, double *ecc,
             double *incl, double *omega, double *argp, double *nu, double *m,
             double *arglat, double *truelon, double *lonper)
{
  double mu = 398600.4418e+9;
  double small = 1e-10;
  double undefined = 999999.1;

  double magr= norm(r);
  double magv= norm(v);

  // ------------------  find h n and e vectors   ----------------
  double hbar[3];
  cross(r, v, hbar);
  double magh = norm(hbar);

  if (magh > small)
  {
    double nbar[3];
    nbar[0] = -hbar[1];
    nbar[1] = hbar[0];
    nbar[2] = 0.0;
    double magn = norm(nbar);
    double c1 = magv*magv - mu/magr;
    double rdotv = dot(r, v);
    double ebar[3];
    for (int i = 0; i < 3; i++)
    {
      ebar[i] = (c1*r[i] - rdotv*v[i])/mu;
    }
    *ecc = norm(ebar);

    // ------------  find a e and semi-latus rectum   ----------
    double sme = (magv*magv*0.5) - (mu/magr);
    *a = -mu/(2.0*sme);
    *p = magh*magh/mu;

    // -----------------  find inclination   -------------------
    double hk= hbar[2]/magh;
    *incl = acos(hk);

    // --------  determine type of orbit for later use  --------
    // ------ elliptical, parabolic, hyperbolic inclined -------
    char typeorbit[2] = {'e', 'i'};
    if (*ecc < small)
    {
      if (*incl < small || fabs(*incl - M_PI) < small)
      {
        // ----------------  circular equatorial ---------------
        typeorbit[0] = 'c';
        typeorbit[1] = 'e';
      }
      else
      {
        // ----------------  circular inclined -----------------
        typeorbit[0] = 'c';
        typeorbit[1] = 'i';
      }
    }
    else
    {
      // ----- elliptical, parabolic, hyperbolic equatorial ----
      if (*incl < small || fabs(*incl - M_PI) < small)
      {
        typeorbit[0] = 'e';
        typeorbit[1] = 'e';
      }
    }

    // ----------  find longitude of ascending node ------------
    double temp;
    if (magn > small)
    {
      temp = nbar[0]/magn;
      if (fabs(temp) > 1.0)
      {
        temp = (double) sign(temp);
      }
      *omega = acos(temp);
      if (nbar[1] < 0.0)
      {
        *omega = 2*M_PI - *omega;
      }
    }
    else
    {
      *omega = undefined;
    }

    // ---------------- find argument of perigee ---------------
    if (typeorbit[0] == 'e' && typeorbit[1] == 'i')
    {
      *argp = angl(nbar, ebar);
      if(ebar[2] < 0.0)
      {
        *argp = 2*M_PI - *argp;
      }
    }
    else
    {
      *argp = undefined;
    }

    // ------------  find true anomaly at epoch    -------------
    if (typeorbit[0] == 'e')
    {
      *nu = angl(ebar, r);
      if(rdotv < 0.0)
      {
        *nu = 2*M_PI - *nu;
      }
    }
    else
    {
      *nu = undefined;
    }

    // ----  find argument of latitude - circular inclined -----
    if (typeorbit[0] == 'c' && typeorbit[1] == 'i')
    {
      *arglat = angl(nbar, r);
      if(r[2] < 0.0)
      {
        *arglat = 2*M_PI - *arglat;
      }
      *m = *arglat;
    }
    else
    {
      *arglat = undefined;
    }

    // -- find longitude of perigee - elliptical equatorial ----
    if (*ecc > small && typeorbit[0] == 'e' && typeorbit[1] == 'e')
    {
      temp = ebar[0] / *ecc;
      if (fabs(temp) > 1.0)
      {
        temp = (double) sign(temp);
      }
      *lonper = acos(temp);
      if (ebar[1] < 0.0)
      {
        *lonper = 2*M_PI - *lonper;
      }
    }
    else
    {
      *lonper = undefined;
    }

    // -------- find true longitude - circular equatorial ------
    if (magr > small && typeorbit[0] == 'c' && typeorbit[1] == 'e')
    {
      temp = r[0]/magr;
      if (fabs(temp) > 1.0)
      {
        temp = (double) sign(temp);
      }
      *truelon = acos(temp);
      if (r[1] < 0.0)
      {
        *truelon = 2*M_PI - *truelon;
      }
      *m = *truelon;
    }
    else
    {
      *truelon = undefined;
    }

    // ------------ find mean anomaly for all orbits -----------
    if (typeorbit[0] == 'e')
    {
      double e;
      newtonnu(*ecc, *nu, &e, m);
    }
  }
  else
  {
    *p       = undefined;
    *a       = undefined;
    *ecc     = undefined;
    *incl    = undefined;
    *omega   = undefined;
    *argp    = undefined;
    *nu      = undefined;
    *m       = undefined;
    *arglat  = undefined;
    *truelon = undefined;
    *lonper  = undefined;
  }
}
//------------------------------------------------------------------------------
