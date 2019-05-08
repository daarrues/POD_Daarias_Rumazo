//------------------------------------------------------------------------------
//                              MatLabUtils
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file MatLabUtils.h
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/04/17
 *
 * Este fichero contiene las cabeceras para las
 * funciones de alto nivel de MatLab necesarias
 * para este proyecto.
 */
//------------------------------------------------------------------------------
#ifndef MATLABUTILS_H
#define MATLABUTILS_H

double norm(double[3]);
double dot(double[3], double[3]);
void cross(double[3], double[3], double[3]);
void zeros(double[], int);
void prodMatr(double[3][3], double[3][3], double[3][3]);
void trans(double[3][3], double[3][3]);
double det(double[3][3]);
int roots(double[], int, double[]);
void unit(double[3], double[3]);
int sign(double);
int fix(double);

#endif // MATLABUTILS_H
