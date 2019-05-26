//------------------------------------------------------------------------------
//                                    newtonnu
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file newtonnu.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero newtonnu.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include <math.h>

//------------------------------------------------------------------------------
//  void newtonnu(double ecc, double nu, double *e0, double *m)
//------------------------------------------------------------------------------
/**
 * newtonnu: Solves keplers equation when the true anomaly is known.
 *   the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
 *   the parabolic limit at 168 is arbitrary. the hyperbolic anomaly is also
 *   limited. the hyperbolic sine is used because it's not double valued.
 *
 * Inputs:
 * @parm <ecc> eccentricity
 * @parm <nu> true anomaly
 *
 * Output:
 * @return <e0> eccentric anomaly
 * @return <m> mean anomaly
 */
//------------------------------------------------------------------------------
void newtonnu(double ecc, double nu, double *e0, double *m)
{
  *e0 = 999999.9;
  *m = 999999.9;
  double small = 0.00000001;
  double sine, cose;

  // --------------------------- circular ------------------------
  if (fabs( ecc ) < small)
  {
    *m = nu;
    *e0 = nu;
  }
  else
  {
    // ---------------------- elliptical -----------------------
    if (ecc < 1.0-small )
    {
      sine = ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
      cose = ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
      *e0 = atan2(sine, cose);
      *m = *e0 - ecc*sin(*e0);
    }
    else
    {
      // -------------------- hyperbolic  --------------------
      if ( ecc > 1.0 + small  )
      {
        if ((ecc > 1.0 ) && (fabs(nu)+0.00001 < M_PI-acos(1.0 /ecc)))
        {
          sine = ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
          *e0 = asinh(sine);
          *m = ecc*sinh(*e0) - (*e0);
        }
      }
      else
      {
        // ----------------- parabolic ---------------------
        if ( fabs(nu) < 168.0*M_PI/180.0  )
        {
          *e0 = tan(nu*0.5);
          *m = *e0 + ((*e0)*(*e0)*(*e0))/3.0;
        }
      }
    }
  }

  if ( ecc < 1.0  )
  {
    *m = fmod(*m,2.0*M_PI);
    if ( *m < 0.0  )
    {
      *m = *m + 2.0 *M_PI;
    }
    *e0 = fmod(*e0,2.0*M_PI);
  }
}
