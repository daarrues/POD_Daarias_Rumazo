//------------------------------------------------------------------------------
//                                   IERS
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file IERS.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/23
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero IERS.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "SAT_Const.h"
#include <math.h>

//------------------------------------------------------------------------------
//  IERS(double *eop[13], double Mjd_UTC, char interp, double *UT1_UTC,
//  double *TAI_UTC, double *x_pole, double *y_pole, double *ddpsi, double *ddeps)
//------------------------------------------------------------------------------
/**
 * Management of IERS time and polar motion data.
 *
 * @param <eop> File data (in).
 * @param <nop> number of rows (in).
 * @param <Mjd_UTC> (in).
 * @param <interp> (in).
 * @param <UT1_UTC> (out).
 * @param <TAI_UTC> (out).
 * @param <x_pole> (out).
 * @param <y_pole> (out).
 * @param <ddpsi> (out).
 * @param <ddeps> (out).
 */
//------------------------------------------------------------------------------
void IERS(double *eop[13], int nop, double Mjd_UTC, char interp, double *UT1_UTC,
  double *TAI_UTC, double *x_pole, double *y_pole, double *ddpsi, double *ddeps)
{
  if(interp == 'l')
  {
    // linear interpolation
    int mj = floor(Mjd_UTC);

    double preeop[13];
    double nexteop[13];
    for(int j = 0; j < nop; j++)
    {
      if(mj == (int) eop[3][j])
      {
        for(int i = 0; i < 13; i++)
        {
          preeop[i] = eop[i][j];
          nexteop[i] = eop[i][j+1];
        }
        break;
      }
    }

    double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
    double fixf = mfme / 1440;

    // Setting of IERS Earth rotation parameters
    // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    *UT1_UTC = preeop[6] + (nexteop[6] - preeop[6])*fixf;
    *TAI_UTC = preeop[12];
    *x_pole  = preeop[4] + (nexteop[4] - preeop[4])*fixf;
    *y_pole  = preeop[5] + (nexteop[5] - preeop[5])*fixf;
    *ddpsi   = preeop[8] + (nexteop[8] - preeop[8])*fixf;
    *ddeps   = preeop[9] + (nexteop[9] - preeop[9])*fixf;

    *x_pole  = *x_pole/Arcs;  // Pole coordinate [rad]
    *y_pole  = *y_pole/Arcs;  // Pole coordinate [rad]
    *ddpsi   = *ddpsi/Arcs;
    *ddeps   = *ddeps/Arcs;
  }
  else if(interp == 'n')
  {
    int mj = floor(Mjd_UTC);

    double myeop[13];
    for(int j = 0; j < nop; j++)
    {
      if(mj == (int) eop[3][j])
      {
        for(int i = 0; i < 13; i++)
        {
          myeop[i] = eop[i][j];
        }
        break;
      }
    }

    // Setting of IERS Earth rotation parameters
    // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    *UT1_UTC = myeop[6];      // UT1-UTC time difference [s]
    *TAI_UTC = myeop[12];     // TAI-UTC time difference [s]
    *x_pole  = myeop[4]/Arcs; // Pole coordinate [rad]
    *y_pole  = myeop[5]/Arcs; // Pole coordinate [rad]
    *ddpsi   = myeop[8]/Arcs;
    *ddeps   = myeop[9]/Arcs;
  }
}
