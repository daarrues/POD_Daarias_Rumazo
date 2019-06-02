//------------------------------------------------------------------------------
//                              TestDaarias
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file testDaarias.h
 * @author Daniel Arias Ruiz-Esquide
 * @date Created: 2019/05/15
 *
 * Este fichero contiene las cabeceras de los
 * test unitarios para las funciones del proyecto
 * implementadas por Daniel Arias Ruiz-Esquide.
 */
//------------------------------------------------------------------------------
#ifndef TEST_DAARIAS_H
#define TEST_DAARIAS_H

typedef int bool;

void testD8rt1(bool);
void testD8rt2(bool);
void testD8rt3(bool);

void testTlamb1(bool);
void testTlamb2(bool);
void testTlamb3(bool);

void testXlamb1(bool);
void testXlamb2(bool);
void testXlamb3(bool);

void testVlamb1(bool);
void testVlamb2(bool);
void testVlamb3(bool);

void testLambert_gooding1(bool);
void testLambert_gooding2(bool);
void testLambert_gooding3(bool);

void testR_x1(bool);
void testR_x2(bool);
void testR_x3(bool);

void testR_y1(bool);
void testR_y2(bool);
void testR_y3(bool);

void testR_z1(bool);
void testR_z2(bool);
void testR_z3(bool);

void testAngl1(bool);
void testAngl2(bool);
void testAngl3(bool);

void testMjday1(bool);
void testMjday2(bool);
void testMjday3(bool);

void testNewtonnu1(bool);
void testNewtonnu2(bool);
void testNewtonnu3(bool);

void testPosition1(bool);
void testPosition2(bool);
void testPosition3(bool);

void testGibbs1(bool);
void testGibbs2(bool);
void testGibbs3(bool);

void testHgibbs1(bool);
void testHgibbs2(bool);
void testHgibbs3(bool);

void testAnglesg(bool);

#endif // TEST_DAARIAS_H
