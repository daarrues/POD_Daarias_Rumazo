//------------------------------------------------------------------------------
//                              TestRumazo
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file testRumazo.h
 * @author Rubén Mazo Tomás
 * @date Created: 2019/05/15
 *
 * Este fichero contiene las cabeceras de los
 * test unitarios para las funciones del proyecto
 * implementadas por Rubén Mazo Tomás.
 */
//------------------------------------------------------------------------------
#ifndef TEST_RUMAZO_H
#define TEST_RUMAZO_H

typedef int bool;

void testMeanObliquity1(bool);
void testMeanObliquity2(bool);
void testMeanObliquity3(bool);

void testNutAngles1(bool);
void testNutAngles2(bool);
void testNutAngles3(bool);

void testEqnEquinox1(bool);
void testEqnEquinox2(bool);
void testEqnEquinox3(bool);

void testFrac1(bool);
void testFrac2(bool);
void testFrac3(bool);

void testGmst1(bool);
void testGmst2(bool);
void testGmst3(bool);

void testTimediff1(bool);
void testTimediff2(bool);
void testTimediff3(bool);

void testIERS1(bool);
void testIERS2(bool);
void testIERS3(bool);

void testGast1(bool);
void testGast2(bool);
void testGast3(bool);

void testGHAMatrix1(bool);
void testGHAMatrix2(bool);
void testGHAMatrix3(bool);

void testPoleMatrix1(bool);
void testPoleMatrix2(bool);
void testPoleMatrix3(bool);

void testNutMatrix1(bool);
void testNutMatrix2(bool);
void testNutMatrix3(bool);

void testPrecMatrix1(bool);
void testPrecMatrix2(bool);
void testPrecMatrix3(bool);

void testDoubler1(bool);
void testDoubler2(bool);
void testDoubler3(bool);

void testRv2coe1(bool);
void testRv2coe2(bool);
void testRv2coe3(bool);

void testAnglesdr1(bool);
void testAnglesdr2(bool);
void testAnglesdr3(bool);

#endif // TEST_RUMAZO_H
