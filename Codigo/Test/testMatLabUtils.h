//------------------------------------------------------------------------------
//                              TestMatLabUtils
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file TestMatLabUtils.h
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/04/17
 *
 * Este fichero contiene las cabeceras de los
 * test unitarios para las funciones de alto
 * nivel de MatLab necesarias para este proyecto.
 */
//------------------------------------------------------------------------------
#ifndef TEST_MATLABUTILS_H
#define TEST_MATLABUTILS_H

typedef int bool;

void testNorm1(bool);
void testNorm2(bool);
void testNorm3(bool);

void testDot1(bool);
void testDot2(bool);
void testDot3(bool);

void testCross1(bool);
void testCross2(bool);
void testCross3(bool);

void testZeros(bool);

void testProdMatr1(bool);
void testProdMatr2(bool);
void testProdMatr3(bool);

void testTrans(bool);

void testDet(bool);

void testRoots(bool);

void testUnit1(bool);
void testUnit2(bool);
void testUnit3(bool);

void testSign1(bool);
void testSign2(bool);
void testSign3(bool);

void testFix1(bool);
void testFix2(bool);
void testFix3(bool);

void testAll1(bool);
void testAll2(bool);
void testAll3(bool);

#endif // TEST_MATLABUTILS_H
