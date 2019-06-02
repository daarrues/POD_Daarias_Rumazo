//------------------------------------------------------------------------------
//                              Test
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file test.h
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/02
 *
 * Este fichero contiene la implementación para la función
 * que realiza todos los test.
 */
//------------------------------------------------------------------------------
#include "test.h"
#include "testMatLabUtils.h"
#include "testDaarias.h"
#include "testRumazo.h"
#include <stdio.h>

typedef int bool;
#define true 1
#define false 0

//------------------------------------------------------------------------------
//  void test()
//------------------------------------------------------------------------------
/**
 *  Función que pasa todos los test unitarios
 */
//------------------------------------------------------------------------------
void test()
{
  printf("TEST MATLABUTILS\n\n");

  // Test norm
  printf("Probando norm!\n");
  testNorm1(false);
  testNorm2(false);
  testNorm3(false);
  printf("norm finalizado!\n\n");

  // Test dot
  printf("Probando dot!\n");
  testDot1(false);
  testDot2(false);
  testDot3(false);
  printf("dot finalizado!\n\n");

  // Test Cross
  printf("Probando cross!\n");
  testCross1(false);
  testCross2(false);
  testCross3(false);
  printf("cross finalizado!\n\n");

  // Test zeros
  printf("Probando zeros!\n");
  testZeros(false);
  printf("zeros finalizado!\n\n");

  // Test prodMatr
  printf("Probando prodMatr!\n");
  testProdMatr1(false);
  testProdMatr2(false);
  testProdMatr3(false);
  printf("prodMatr finalizado!\n\n");

  // Test trans
  printf("Probando trans!\n");
  testTrans(false);
  printf("trans finalizado!\n\n");

  // Test det
  printf("Probando det!\n");
  testDet(false);
  printf("det finalizado!\n\n");

  // Test roots
  printf("Probando roots!\n");
  testRoots(false);
  printf("roots finalizado!\n\n");

  // Test unit
  printf("Probando unit!\n");
  testUnit1(false);
  testUnit2(false);
  testUnit3(false);
  printf("unit finalizado!\n\n");

  // Test sign
  printf("Probando sign!\n");
  testSign1(false);
  testSign2(false);
  testSign3(false);
  printf("sign finalizado!\n\n");

  // Test fix
  printf("Probando fix!\n");
  testFix1(false);
  testFix2(false);
  testFix3(false);
  printf("fix finalizado!\n\n");

  // Test all
  printf("Probando all!\n");
  testAll1(false);
  testAll2(false);
  testAll3(false);
  printf("all finalizado!\n\n");

  printf("TEST PRELIMINARY_ORBIT_DETERMINATION\n\n");

  // Test MeanObliquity
	printf("Probando MeanObliquity!\n");
	testMeanObliquity1(false);
	testMeanObliquity2(false);
	testMeanObliquity3(false);
	printf("MeanObliquity finalizado!\n\n");

	// Test NutAngles
	printf("Probando NutAngles!\n");
	testNutAngles1(false);
	testNutAngles2(false);
	testNutAngles3(false);
	printf("NutAngles finalizado!\n\n");

	// Test EqnEquinox
	printf("Probando EqnEquinox!\n");
	testEqnEquinox1(false);
	testEqnEquinox2(false);
	testEqnEquinox3(false);
	printf("EqnEquinox finalizado!\n\n");

	// Test Frac
	printf("Probando Frac!\n");
	testFrac1(false);
	testFrac2(false);
	testFrac3(false);
	printf("Frac finalizado!\n\n");

	// Test gmst
	printf("Probando gmst!\n");
	testGmst1(false);
	testGmst2(false);
	testGmst3(false);
	printf("gmst finalizado!\n\n");

	// Test timediff
	printf("Probando timediff!\n");
	testTimediff1(false);
	testTimediff2(false);
	testTimediff3(false);
	printf("timediff finalizado!\n\n");

	// Test IERS
	printf("Probando IERS!\n");
	testIERS1(false);
	testIERS2(false);
	testIERS3(false);
	printf("IERS finalizado!\n\n");

	// Test gast
	printf("Probando gast!\n");
	testGast1(false);
	testGast2(false);
	testGast3(false);
	printf("gast finalizado!\n\n");

  // Test R_x
	printf("Probando R_x!\n");
	testR_x1(false);
	testR_x2(false);
	testR_x3(false);
	printf("R_x finalizado!\n\n");

	// Test R_y
	printf("Probando R_y!\n");
	testR_y1(false);
	testR_y2(false);
	testR_y3(false);
	printf("R_y finalizado!\n\n");

	// Test R_z
	printf("Probando R_z!\n");
	testR_z1(false);
	testR_z2(false);
	testR_z3(false);
	printf("R_z finalizado!\n\n");

	// Test GHAMatrix
	printf("Probando GHAMatrix!\n");
	testGHAMatrix1(false);
	testGHAMatrix2(false);
	testGHAMatrix3(false);
	printf("GHAMatrix finalizado!\n\n");

	// Test PoleMatrix
	printf("Probando PoleMatrix!\n");
	testPoleMatrix1(false);
	testPoleMatrix2(false);
	testPoleMatrix3(false);
	printf("PoleMatrix finalizado!\n\n");

	// Test NutMatrix
	printf("Probando NutMatrix!\n");
	testNutMatrix1(false);
	testNutMatrix2(false);
	testNutMatrix3(false);
	printf("NutMatrix finalizado!\n\n");

	// Test PrecMatrix
	printf("Probando PrecMatrix!\n");
	testPrecMatrix1(false);
	testPrecMatrix2(false);
	testPrecMatrix3(false);
	printf("PrecMatrix finalizado!\n\n");

  // Test Mjday
	printf("Probando Mjday!\n");
	testMjday1(false);
	testMjday2(false);
	testMjday3(false);
	printf("Mjday finalizado!\n\n");

	// Test Position
	printf("Probando Position!\n");
	testPosition1(false);
	testPosition2(false);
	testPosition3(false);
	printf("Position finalizado!\n\n");

  // Test d8rt
	printf("Probando d8rt!\n");
	testD8rt1(false);
	testD8rt2(false);
	testD8rt3(false);
	printf("d8rt finalizado!\n\n");

  // Test tlamb
	printf("Probando tlamb!\n");
	testTlamb1(false);
	testTlamb2(false);
	testTlamb3(false);
	printf("tlamb finalizado!\n\n");

	// Test xlamb
	printf("Probando xlamb!\n");
	testXlamb1(false);
	testXlamb2(false);
	testXlamb3(false);
	printf("xlamb finalizado!\n\n");

	// Test vlamb
	printf("Probando vlamb!\n");
	testVlamb1(false);
	testVlamb2(false);
	testVlamb3(false);
	printf("vlamb finalizado!\n\n");

	// Test lambert_gooding
	printf("Probando lambert_gooding!\n");
	testLambert_gooding1(false);
	testLambert_gooding2(false);
	testLambert_gooding3(false);
	printf("lambert_gooding finalizado!\n\n");

	// Test doubler
	printf("Probando doubler!\n");
	testDoubler1(false);
	testDoubler2(false);
	testDoubler3(false);
	printf("doubler finalizado!\n\n");

  // Test angl
	printf("Probando angl!\n");
	testAngl1(false);
	testAngl2(false);
	testAngl3(false);
	printf("angl finalizado!\n\n");

  // Test newtonnu
	printf("Probando newtonnu!\n");
	testNewtonnu1(false);
	testNewtonnu2(false);
	testNewtonnu3(false);
	printf("newtonnu finalizado!\n\n");

	// Test rv2coe
	printf("Probando rv2coe!\n");
	testRv2coe1(false);
	testRv2coe2(false);
	testRv2coe3(false);
	printf("rv2coe finalizado!\n\n");

  // Test gibbs
	printf("Probando gibbs!\n");
	testGibbs1(false);
	testGibbs2(false);
	testGibbs3(false);
	printf("gibbs finalizado!\n\n");

	// Test hgibbs
	printf("Probando hgibbs!\n");
	testHgibbs1(false);
	testHgibbs2(false);
	testHgibbs3(false);
	printf("hgibbs finalizado!\n\n");

	// Test anglesdr
	printf("Probando anglesdr!\n");
	testAnglesdr1(false);
	testAnglesdr2(false);
	testAnglesdr3(false);
	printf("anglesdr finalizado!\n\n");

  // Test anglesg
	printf("Probando anglesg!\n");
	testAnglesg(false);
	printf("anglesg finalizado!\n\n");

  // Final
	printf("TODOS LOS TEST SUPERADOS\n\n");
}
//------------------------------------------------------------------------------
