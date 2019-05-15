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

void lambert_gooding(double[], double[], double, double, double, double[],
                     double[]);
void vlamb(double, double, double, double, double, double, double[], double[],
           double[], double[]);
void tlamb(double, double, double, double, double, double, double, double,
           double);
double d8rt(double);
void xlamb(double, double, double, double, double, double, double);

#endif // LAMBERT_GOODING_H
