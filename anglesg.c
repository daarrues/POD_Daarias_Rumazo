//------------------------------------------------------------------------------
//                              Anglesg
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file anglesg.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/01
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero anglesg.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "anglesg.h"
#include "MatLabUtils/MatLabUtils.h"
#include <math.h>

//------------------------------------------------------------------------------
//  void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1,
//               double Delta2, double Delta13, double JD1, double JD2,
//               double JD3, double RS1[], double RS2[], double RS3[],
//               double R2[], double V2[])
//------------------------------------------------------------------------------
/**
 *  anglesgauss: Solves the problem of orbit determination using three
 *               optical sightings.
 *
 * Inputs:
 * @param <Alpha1> - right ascension #1 (rad)
 * @param <Alpha2> - right ascension #2 (rad)
 * @param <Alpha3> - right ascension #3 (rad)
 * @param <Delta1> - declination #1 (rad)
 * @param <Delta2> - declination #2 (rad)
 * @param <Delta3> - declination #3 (rad)
 * @param <JD1> - julian date of 1st sighting (days since 4713 bc)
 * @param <JD2> - julian date of 2nd sighting (days since 4713 bc)
 * @param <JD3> - julian date of 3rd sighting (days since 4713 bc)
 * @param <RS1> - ijk site1 position vector (m)
 * @param <RS2> - ijk site2 position vector (m)
 * @param <RS3> - ijk site3 position vector (m)
 *
 *  Outputs:
 * @return <R2> - ijk position vector at t2 (m)
 * @return <V2> - ijk velocity vector at t2 (m/s)
 */
//------------------------------------------------------------------------------
void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1,
             double Delta2, double Delta13, double JD1, double JD2, double JD3,
             double RS1[], double RS2[], double RS3[], double R2[], double V2[])
{
  double Mu = 398600.4418e9;
  double Rad = 180/M_PI;
  double R1[3], R2[3], R3[3];
  zeros(R1,3);
  zeros(R2,3);
  zeros(R3,3);

  //---------- set middle to 0, find deltas to others -----------
  double Tau1 = (JD1-JD2)*86400;
  double Tau3 = (JD3-JD2)*86400;

  //----------------  Find Line of SIGHT vectors  ---------------
  double L1[3], L2[3], L3[3];
  L1[0] = cos(Delta1)*cos(Alpha1);
  L1[1] = cos(Delta1)*sin(Alpha1);
  L1[2] = sin(Delta1);

  L2[0] = cos(Delta2)*cos(Alpha2);
  L2[1] = cos(Delta2)*sin(Alpha2);
  L2[2] = sin(Delta2);

  L3[0] = cos(Delta3)*cos(Alpha3);
  L3[1] = cos(Delta3)*sin(Alpha3);
  L3[2] = sin(Delta3);

  //------------- Find L matrix and determinant -----------------
  //--------- Called LMatI since it is only used for determ -----
  double LMatIi[3][3], RSMat[3][3];
  for (i = 0; i < 3; i++)
  {
    LMatIi[i][0] = L1[i];
    LMatIi[i][1] = L2[i];
    LMatIi[i][2] = L3[i];
    RSMat[i][0] = RS1[i];
    RSMat[i][1] = RS2[i];
    RSMat[i][2] = RS3[i];
  }

  double D = det(LMatIi);

  //------------------ Now assign the inverse -------------------
  double LIR[3][3];
  LMatI[0][0] = ( L2[1]*L3[2]-L2[2]*L3[1])/D;
  LMatI[1][0] = (-L1[1]*L3[2]+L1[2]*L3[1])/D;
  LMatI[2][0] = ( L1[1]*L2[2]-L1[2]*L2[1])/D;
  LMatI[0][1] = (-L2[0]*L3[2]+L2[2]*L3[0])/D;
  LMatI[1][1] = ( L1[0]*L3[2]-L1[2]*L3[0])/D;
  LMatI[2][1] = (-L1[0]*L2[2]+L1[2]*L2[0])/D;
  LMatI[0][2] = ( L2[0]*L3[1]-L2[1]*L3[0])/D;
  LMatI[1][2] = (-L1[0]*L3[1]+L1[1]*L3[0])/D;
  LMatI[2][2] = ( L1[0]*L2[1]-L1[1]*L2[0])/D;

  prodMatr(LMatI, RSMat, LIR);

  //------------ Find f and g series at 1st and 3rd obs ---------
  //      speed by assuming circ sat vel for uDot here ??
  //      some similartities in 1/6t3t1 ...
  //--- keep separated this time ----
  double a1  = Tau3/(Tau3 - Tau1);
  double a1u = (Tau3*((Tau3-Tau1)*(Tau3-Tau1) - Tau3*Tau3 ))/(6.0*(Tau3 - Tau1));
  double a3  = -Tau1 / (Tau3 - Tau1);
  double a3u = -(Tau1*((Tau3-Tau1)*(Tau3-Tau1) - Tau1*Tau1 ))/(6.0*(Tau3 - Tau1));

  //--- Form initial guess of r2 ----
  double D1 = LIR[1][0]*a1 - LIR[1][1] + LIR[1][2]*a3;
  double D2 = LIR[1][0]*a1u + LIR[1][2]*a3u;

  //------- Solve eighth order poly NOT same as LAPLACE ---------
  double L2DotRS= dot(L2,RS2);
  double magRS2 = norm(RS2);
  double Poly[POL_DEG + 1], rootarr[POL_DEG];
  zeros(Poly, POL_DEG + 1);
  Poly[0] = 1.0;  // r2^8th variable!!!!!!!!!!!!!!
  Poly[2] = -(D1*D1 + 2.0*D1*L2DotRS + magRS2^2);
  Poly[5] =  -2.0*Mu*(L2DotRS*D2 + D1*D2);
  Poly[8] =  -Mu*Mu*D2*D2;

  int nRoots = roots(Poly, POL_DEG, rootarr);

  //------------------ Select the correct root ------------------
  double BigR2 = 0.0;

  for (int j = 0; j < nRoots; j++)
  {
    if (rootarr[j] > BigR2)
    {
      BigR2 = rootarr[j];
    }
  }

  //------------ Solve matrix with u2 better known --------------
  double u = Mu/(BigR2*BigR2*BigR2);

  double c1 = a1+a1u*u;
  double c3 = a3+a3u*u;
  double CMat[3], RhoMat[3];
  CMat[0] = -c1;
  CMat[1] = 1.0;
  CMat[2] = -c3;
  for(int i = 0; i < 3; i++)
  {
    dot(LIR[i], CMat, RhoMat[i]);
  }

  // Rhoold1 = RhoMat[0]/c1;
  Rhoold2 = -RhoMat[1];
  // Rhoold3 = RhoMat[2]/c3;

  double ll, Rho2, theta, theta1, Theta, Theta1, copa, p, a, ecc, incl, omega,
  argp, Nu, m, u, l, ArgPer, magR1, magR2, magR3, U, RDot, UDot, TauSqr, f1, f3,
  g1, g3;
  double v1[3], v2[3];
  double R1t[3][3]
  char error[12];
  while ( (fabs(Rhoold2-Rho2)>1e-12) && (ll<=2) )
  {
    ll = ll + 1;
    Rho2 = Rhoold2; // reset now that inside while loop
    //---------- Now form the three position vectors ----------
    for (int i = 0; i < 3; i++)
    {
      R1[i] = RhoMat[0]*L1[i]/c1 + RS1[i];
      R2[i] = -RhoMat[1]*L2[i] + RS2[i];
      R3[i] = RhoMat[2]*L3[i]/c3 + RS3[i];
    }

    gibbs(R1, R2, R3, v2, &theta, &theta1, &copa, error);

    if (!strcmp(error,"          ok") && (copa < 1/Rad))
    {
      //--- HGibbs to get middle vector ----
      hgibbs(R1, R2, R3, JD1, JD2, JD3, v2, &theta, &theta1, &copa, error);
    }
    lambert_gooding(R1, R2, (JD2-JD1)*86400, Mu, 0.0, 1, v1,V2);
    rv2coe(R2,V2, &p, &a, &ecc, &incl, &omega, &argp, &Nu, &m, &u, &l, &ArgPer);
    magR2 = norm(R2);

    if ( ll <= 2 )
    {
      //--- Now get an improved estimate of the f and g series --
      //       .or. can the analytic functions be found now??
      U = Mu/(magR2^3);
      RDot = dot(R2,V2)/magR2;
      UDot = (-3.0*Mu*RDot)/(magR2^4);

      TauSqr= Tau1*Tau1;
      f1 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau1
           + (1.0/24.0) * U*U*TauSqr*TauSqr
           + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau1;
      g1 = Tau1 - (1.0/6.0)*U*Tau1*TauSqr - (1.0/12.0) *
           UDot*TauSqr*TauSqr
           + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau1
           + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;

      TauSqr = Tau3*Tau3;
      f3 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau3
           + (1.0/24.0) * U*U*TauSqr*TauSqr
           + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau3;
      g3 = Tau3 - (1.0/6.0)*U*Tau3*TauSqr - (1.0/12.0) *
           UDot*TauSqr*TauSqr
           + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau3
           + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
    }
    else
    {
      //-------- Now use exact method to find f and g -----------
      Theta = angl(R1,R2);
      Theta1 = angl(R2,R3);
      magR1 = norm(R1);
      magR3 = norm(R3);

      f1 = 1.0 - ( (magR1*(1.0 - cos(Theta))/p ) );
      g1 = ( magR1*magR2*sin(-theta) ) / sqrt(p);  // -ANGLE because backwards!!
      f3 = 1.0 - ( (magR3*(1.0 - cos(Theta1))/p ) );
      g3 = ( magR3*magR2*sin(theta1) )/sqrt(p);
    }

    c1 =  g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);

    //----- Solve for all three ranges via matrix equation ----
    CMat[0] = -c1;
    CMat[1] = 1.0;
    CMat[2] = -c3;
    for(int i = 0; i < 3; i++)
    {
      dot(LIR[i], CMat, RhoMat[i]);
    }

    // Rhoold1 = RhoMat[0]/c1;
    Rhoold2 = -RhoMat[1];
    // Rhoold3 = RhoMat[2]/c3;
    //----------------- Check for convergence -----------------
  } // DO WHILE the ranges are still changing
  //---------------- Find all three vectors ri ------------------
  for (i = 0; i < 3; i++)
  {
    R1[i] = RhoMat[0]*L1[i]/c1 + RS1[i];
    R2[i] = -RhoMat[1]*L2[i] + RS2[i];
    R3[i] = RhoMat[2]*L3[i]/c3 + RS3[i];
  }
}
