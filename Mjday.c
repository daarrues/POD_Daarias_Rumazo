//------------------------------------------------------------------------------
//                                    Mjday
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file Mjday.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/25
 *
 * Este fichero contiene las implementaciones para las
 * funciones del fichero Mjday.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#include "MatLabUtils/MatLabUtils.h"
#include <math.h>
#include <stdio.h>

//------------------------------------------------------------------------------
//  void Mjday(double year, double month, double day, double hour, double min,
//             double sec, double* Mjd)
//------------------------------------------------------------------------------
/**
 * Mjday: Modified Julian Date from calendar date and time
 *
 * Inputs:
 * @parm <year>
 * @parm <month>
 * @parm <day>
 * @parm <hour>
 * @parm <min>
 * @parm <sec>
 *
 * Output:
 * @return <Mjd> Modified Julian Date
 */
//------------------------------------------------------------------------------
void Mjday(double year, double month, double day, double hour, double min,
            double sec, double* Mjd)
{
  double y = year;
  double m = month;
  double a;
  double b = 0.0;
  double c = 0.0;

  if (m <= 2)
  {
    y = y - 1;
    m = m + 12;
  }

  if (y < 0)
  {
    c = -.75;
  }

  // check for valid calendar date
  if (year < 1582)
  {
    // null
  }
  else if (year > 1582)
  {
    a = fix(y / 100);
    b = 2 - a + floor(a / 4);
  }
  else if (month < 10)
  {
    // null
  }
  else if (month > 10)
  {
    a = fix(y / 100);
    b = 2 - a + floor(a / 4);
  }
  else if (day <= 4)
  {
    // null
  }
  else if (day > 14)
  {
    a = fix(y / 100);
    b = 2 - a + floor(a / 4);
  }
  else
  {
    printf("\n\n  this is an invalid calendar date!!\n");
    Mjd = NULL;
  }

  double jd = fix(365.25 * y + c) + fix(30.6001 * (m + 1));
  jd = jd + day + b + 1720994.5;
  jd = jd + (hour+min/60+sec/3600)/24;
  *Mjd = jd - 2400000.5;
}
//------------------------------------------------------------------------------
