//------------------------------------------------------------------------------
//                                  timediff
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file timediff.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/22
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero timediff.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "timediff.h"

//------------------------------------------------------------------------------
//  void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI,
//  double *UTC_GPS, double *UT1_GPS, double *TT_UTC, double *GPS_UTC)
//------------------------------------------------------------------------------
/**
 * Time differences [s]
 *
 * @param <UT1_UTC> UT1-UTC time difference [s] (in).
 * @param <TAI_UTC> TAI-UTC time difference [s] (in).
 * @param <UT1_TAI> UT1-TAI time difference [s] (out).
 * @param <UTC_GPS> UTC-GPS time difference [s] (out).
 * @param <UT1_GPS> UT1-GPS time difference [s] (out).
 * @param <TT_UTC> TT-UTC time difference [s] (out).
 * @param <GPS_UTC> GPS-UTC time difference [s] (out).
 */
//------------------------------------------------------------------------------
void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS,
              double *UT1_GPS, double *TT_UTC, double *GPS_UTC)
{
  double TT_TAI  = 32.184;         //  TT-TAI time difference [s]
  double GPS_TAI = -19.0;          // GPS-TAI time difference [s]
  *UT1_TAI = UT1_UTC - TAI_UTC;    // UT1-TAI time difference [s]
  double UTC_TAI = -TAI_UTC;       // UTC-TAI time difference [s]
  *UTC_GPS = UTC_TAI - GPS_TAI;    // UTC-GPS time difference [s]
  *UT1_GPS = *UT1_TAI - GPS_TAI;   // UT1-GPS time difference [s]
  *TT_UTC  = TT_TAI - UTC_TAI;     //  TT-UTC time difference [s]
  *GPS_UTC = GPS_TAI - UTC_TAI;    // GPS-UTC time difference [s]
}
//------------------------------------------------------------------------------
