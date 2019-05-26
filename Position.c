//------------------------------------------------------------------------------
//                              Position
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file Position.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero Position.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "SAT_Const.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void Postion(double lon, double lat, double h, double r[])
//------------------------------------------------------------------------------
/*
 * Position: Position vector (r [m]) from geodetic coordinates
 *           (Longitude [rad], latitude [rad], altitude [m])
 *
 * @param <lon> longitude
 * @param <lat> latitude
 * @param <h> altitude
 *
 * @return <r> position vector
 */
//------------------------------------------------------------------------------
void Postion(double lon, double lat, double h, double r[])
{
  double R_equ = R_Earth;
  double f     = f_Earth;

  double e2     = f*(2-f);  // Square of eccentricity
  double CosLat = cos(lat); // Cosine of geodetic latitude
  double SinLat = sin(lat); // Sine of geodetic latitude

  // Position vector
  double N = R_equ/sqrt(1-e2*SinLat*SinLat);

  r[0] =  (       N+h)*CosLat*cos(lon);
  r[1] =  (       N+h)*CosLat*sin(lon);
  r[2] =  ((1-e2)*N+h)*SinLat;
}
//------------------------------------------------------------------------------
