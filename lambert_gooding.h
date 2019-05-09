//------------------------------------------------------------------------------
//                              lambert_gooding
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
* @file lambert_gooding.h
* @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
* @date Created: 2019/05/09
 *
 * Este fichero contiene las cabeceras para las
 * funciones del fichero lambert_gooding.m (M. Mahooti)
 */
//------------------------------------------------------------------------------

#ifndef LAMBERT_GOODING_H
#define LAMBERT_GOODING_H

void lambert_gooding(double, double, double, int, int, double[][], double[][]);
void vlamb(double, double, double, double, double, int, double[], double[],
           double[], double[]);
void tlamb(m,q,qsqfm1,x,n)
double d8rt(double);
xlamb(m,q,qsqfm1,tin);

#endif // LAMBERT_GOODING_H
