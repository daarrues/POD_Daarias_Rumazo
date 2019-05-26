//------------------------------------------------------------------------------
//                                    angl
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file angl.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero angl.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "MatLabUtils/MatLabUtils.h"
#include <math.h>

//------------------------------------------------------------------------------
//  double angl(double vec1[], double vec2[])
//------------------------------------------------------------------------------
/**
 * Calculates the angle between two vectors [-pi,pi]
 *
 * Inputs:
 * @parm <vec1> vector 1
 * @parm <vec2> vector 2
 *
 * Output:
 * @return <theta> angle between the two vectors [-pi,pi]
 */
//------------------------------------------------------------------------------
double angl(double vec1[], double vec2[])
{
  double small     = 0.00000001;
  double undefined = 999999.1;

  double magv1 = norm(vec1);
  double magv2 = norm(vec2);

  if (magv1*magv2 > small^2)
  {
    double temp= dot(vec1,vec2)/(magv1*magv2);
    if (fabs(temp) > 1)
    {
      temp = sign(temp);
    }
    return acos(temp);
  }
  else
  {
    return undefined;
  }
}
//------------------------------------------------------------------------------
