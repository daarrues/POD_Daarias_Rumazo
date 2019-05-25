//------------------------------------------------------------------------------
//                                   IERS
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file IERS.h
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/05/23
 *
 * Este fichero contiene las cabeceras para las
 * funciones del fichero IERS.m (M. Mahooti)
 */
//------------------------------------------------------------------------------
#ifndef IERS_H
#define IERS_H

void IERS(double*[13], int, double, char,
          double*, double*, double*, double*, double*, double*);

#endif // IERS_H
