//------------------------------------------------------------------------------
//                                 Example6
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file example6.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/02
 *
 * Este fichero contiene la implementación de example6.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "../MatLabUtils/MatLabUtils.h"
#include "anglesdr.h"
#include "examples.h"
#include "EOPDATA.h"
#include "GHAMatrix.h"
#include "IERS.h"
#include "Mjday.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "Position.h"
#include "PrecMatrix.h"
#include "SAT_Const.h"
#include "timediff.h"
#include <stdio.h>

//------------------------------------------------------------------------------
//  void example6()
//------------------------------------------------------------------------------
/**
 *  Preliminary Orbit Determination using Double-R-Iteration method
 *
 * @note References:
 *   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
 *   Applications", Springer Verlag, Heidelberg, 2005.
 *
 *   D. Vallado, "Fundamentals of Astrodynamics and Applications",
 *   3rd Edition, 2007.
 *
 *   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
 */
//------------------------------------------------------------------------------
void example6()
{
  double *eopdata[13];
  leerFichero(eopdata);

  double Y, M, D, h, m, s, rtasc, decl;
  double obs[3][3];

  FILE *fich = fopen("PreliminaryOrbitDetermination/sat6.txt", "r");

  int i = 0;
	while(fscanf(fich, "%lf/%lf/%lf %lf:%lf:%lf %lf %lf\n",
				&Y, &M, &D, &h, &m, &s, &rtasc, &decl) != EOF)
	{
    obs[i][0] = Mjday(Y,M,D,h,m,s);
    obs[i][1] = Rad*rtasc;
    obs[i][2] = Rad*decl;
		i++;
	}

	fclose(fich);

  // station
  double lat = Rad*30.5724;     // (rad)
  double lon = Rad*(-86.2143);  // (rad)
  double alt = 0.0;             // (m)

  double Rs[3];
  Position(lon, lat, alt,Rs);

  double Mjd1 = obs[0][0];
  double Mjd2 = obs[1][0];
  double Mjd3 = obs[2][0];
  double Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole, ddpsi, ddeps, UT1_TAI, UTC_GPS,
         UT1_GPS, TT_UTC, GPS_UTC, Mjd_TT, Mjd_UT1;
  // Original matrix to use
  double P[3][3], N[3][3], E[3][3];
  // Added while migrating to C
  double poleM[3][3], ghaM[3][3], pgM[3][3], pgnM[3][3], Et[3][3];

  // First site
  Mjd_UTC = Mjd1;
  IERS(eopdata, 20026, Mjd_UTC,'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole,
       &ddpsi, &ddeps);
  timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
  Mjd_TT = Mjd_UTC + TT_UTC/86400;
  Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

  PrecMatrix(MJD_J2000, Mjd_TT, P);
  NutMatrix(Mjd_TT, N);
  PoleMatrix(x_pole, y_pole, poleM);
  GHAMatrix(Mjd_UT1, ghaM);
  prodMatr(poleM, ghaM, pgM);
  prodMatr(pgM, N, pgnM);
  prodMatr(pgnM, P, E);
  trans(E, Et);
  double rsite1[3];
  for(int i = 0; i < 3; i++)
  {
    rsite1[i] = dot(Et[i], Rs);
  }

  // Second site
  Mjd_UTC = Mjd2;
  IERS(eopdata, 20026, Mjd_UTC,'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole,
       &ddpsi, &ddeps);
  timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
  Mjd_TT = Mjd_UTC + TT_UTC/86400;
  Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

  PrecMatrix(MJD_J2000, Mjd_TT, P);
  NutMatrix(Mjd_TT, N);
  PoleMatrix(x_pole, y_pole, poleM);
  GHAMatrix(Mjd_UT1, ghaM);
  prodMatr(poleM, ghaM, pgM);
  prodMatr(pgM, N, pgnM);
  prodMatr(pgnM, P, E);
  trans(E, Et);
  double rsite2[3];
  for(int i = 0; i < 3; i++)
  {
    rsite2[i] = dot(Et[i], Rs);
  }

  // Third site
  Mjd_UTC = Mjd3;
  IERS(eopdata, 20026, Mjd_UTC,'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole,
       &ddpsi, &ddeps);
  timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
  Mjd_TT = Mjd_UTC + TT_UTC/86400;
  Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

  PrecMatrix(MJD_J2000, Mjd_TT, P);
  NutMatrix(Mjd_TT, N);
  PoleMatrix(x_pole, y_pole, poleM);
  GHAMatrix(Mjd_UT1, ghaM);
  prodMatr(poleM, ghaM, pgM);
  prodMatr(pgM, N, pgnM);
  prodMatr(pgnM, P, E);
  trans(E, Et);
  double rsite3[3];
  for(int i = 0; i < 3; i++)
  {
    rsite3[i] = dot(Et[i], Rs);
  }

  double r2[3], v2[3];
  anglesdr(obs[0][1], obs[1][1], obs[2][1], obs[0][2], obs[1][2], obs[2][2],
          Mjd1, Mjd2, Mjd3, rsite1, rsite2, rsite3, r2, v2);
  printf("\nDouble-R-Iteration method");
  printf("\nY_apr=\n\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n",
          r2[0]*1e-3, r2[1]*1e-3, r2[2]*1e-3,
          v2[0]*1e-3, v2[1]*1e-3, v2[2]*1e-3);
}
