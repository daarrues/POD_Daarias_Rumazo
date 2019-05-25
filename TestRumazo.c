#include "EOPDATA.h"

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "EqnEquinox.h"
#include "Frac.h"
#include "gmst.h"
#include "timediff.h"
#include "IERS.h"
#include "gast.h"

typedef int bool;
#define true 1
#define false 0

//------------------------------------------------------------------------------
//  void testMeanObliquity1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función MeanObliquity
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testMeanObliquity1(bool verbose)
{
	double n = 54977.6676696643;
	double esperado = 0.409071470558628;
	double obtenido = MeanObliquity(n);

	if(verbose)
	{
		printf("MeanObliquity1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("MeanObliquity1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testMeanObliquity2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función MeanObliquity
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testMeanObliquity2(bool verbose)
{
	double n = 55565.5429708796;
	double esperado = 0.409067817509821;
	double obtenido = MeanObliquity(n);

	if(verbose)
	{
		printf("MeanObliquity2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("MeanObliquity2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testMeanObliquity3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función MeanObliquity
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testMeanObliquity3(bool verbose)
{
	double n = 55565.9051733796;
	double esperado = 0.409067815259099;
	double obtenido = MeanObliquity(n);

	if(verbose)
	{
		printf("MeanObliquity3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("MeanObliquity3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutAngles1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función NutAngles
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNutAngles1(bool verbose)
{
	double n = 54977.6676696643;
	double esperado1 = 6.48693387393849e-05;
	double esperado2 = 2.23051336707266e-05;
	double obtenido1;
	double obtenido2;

	NutAngles(n, &obtenido1, &obtenido2);

	if(verbose)
	{
		printf("NutAngles1:\n");
		printf("Esperado 1: %.20lf \n", esperado1);
		printf("Obtenido 1: %.20lf \n", obtenido1);
		printf("Esperado 2: %.20lf \n", esperado2);
		printf("Obtenido 2: %.20lf \n", obtenido2);
	}

	assert(fabs(esperado1 - obtenido1) < 10e-12);
	assert(fabs(esperado2 - obtenido2) < 10e-12);
	printf("NutAngles1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutAngles2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función NutAngles
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNutAngles2(bool verbose)
{
	double n = 53989.1991812154;
	double esperado1 = 7.09833624503762e-06;
	double esperado2 = 4.67343559525024e-05;
	double obtenido1;
	double obtenido2;

	NutAngles(n, &obtenido1, &obtenido2);

	if(verbose)
	{
		printf("NutAngles2:\n");
		printf("Esperado 1: %.20lf \n", esperado1);
		printf("Obtenido 1: %.20lf \n", obtenido1);
		printf("Esperado 2: %.20lf \n", esperado2);
		printf("Obtenido 2: %.20lf \n", obtenido2);
	}

	assert(fabs(esperado1 - obtenido1) < 10e-12);
	assert(fabs(esperado2 - obtenido2) < 10e-12);
	printf("NutAngles2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutAngles3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función NutAngles
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNutAngles3(bool verbose)
{
	double n = 54332.4868655555;
	double esperado1 = 3.52072380612286e-05;
	double esperado2 = 4.18260166677631e-05;
	double obtenido1;
	double obtenido2;

	NutAngles(n, &obtenido1, &obtenido2);

	if(verbose)
	{
		printf("NutAngles3:\n");
		printf("Esperado 1: %.20lf \n", esperado1);
		printf("Obtenido 1: %.20lf \n", obtenido1);
		printf("Esperado 2: %.20lf \n", esperado2);
		printf("Obtenido 2: %.20lf \n", obtenido2);
	}

	assert(fabs(esperado1 - obtenido1) < 10e-12);
	assert(fabs(esperado2 - obtenido2) < 10e-12);
	printf("NutAngles3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testEqnEquinox1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función EqnEquinox
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testEqnEquinox1(bool verbose)
{
	double n = 54977.6676696643;
	double esperado = 5.95170051422054e-05;
	double obtenido = EqnEquinox(n);

	if(verbose)
	{
		printf("EqnEquinox1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("EqnEquinox1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testEqnEquinox2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función EqnEquinox
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testEqnEquinox2(bool verbose)
{
	double n = 55565.9051733796;
	double esperado = 8.02104092023363e-05;
	double obtenido = EqnEquinox(n);

	if(verbose)
	{
		printf("EqnEquinox2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("EqnEquinox2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testEqnEquinox3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función EqnEquinox
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testEqnEquinox3(bool verbose)
{
	double n = 53989.1991812154;
	double esperado = 6.51263906819133e-06;
	double obtenido = EqnEquinox(n);

	if(verbose)
	{
		printf("EqnEquinox3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("EqnEquinox3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testFrac1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función Frac
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testFrac1(bool verbose)
{
	double n = 10.3456636914056;
	double esperado = 0.3456636914056;
	double obtenido = Frac(n);

	if(verbose)
	{
		printf("Frac1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("Frac1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testFrac2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función Frac
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testFrac2(bool verbose)
{
	double n = 7.17084861050766;
	double esperado = 0.17084861050766;
	double obtenido = Frac(n);

	if(verbose)
	{
		printf("Frac2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("Frac2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testFrac3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función Frac
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testFrac3(bool verbose)
{
	double n = 11.2920772238043;
	double esperado = 0.2920772238043;
	double obtenido = Frac(n);

	if(verbose)
	{
		printf("Frac3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("Frac3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGmst1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función gmst
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGmst1(bool verbose)
{
	double n = 54977.6669066321;
	double esperado = 2.17186902706532;
	double obtenido = gmst(n);

	if(verbose)
	{
		printf("gmst1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-8);
	printf("gmst1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGmst2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función gmst
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGmst2(bool verbose)
{
	double n = 55565.9044057253;
	double esperado = 1.21707647675469;
	double obtenido = gmst(n);

	if(verbose)
	{
		printf("gmst2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-8);
	printf("gmst2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGmst3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función gmst
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGmst3(bool verbose)
{
	double n = 54332.4861092138;
	double esperado = 2.50334500393596;
	double obtenido = gmst(n);

	if(verbose)
	{
		printf("gmst3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-8);
	printf("gmst3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testTimediff1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función timediff
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testTimediff1(bool verbose)
{
	double n1 = 0.258022690875596;
	double n2 = 34.0;
	double esp1 = -33.7419773091244;
	double esp2 = -15.0;
	double esp3 = -14.7419773091244;
	double esp4 = 66.184;
	double esp5 = 15.0;
	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;

	timediff(n1, n2, &obt1, &obt2, &obt3, &obt4, &obt5);

	if(verbose)
	{
		printf("timediff1:\n");
		printf("Esperado 1: %.20lf \n", esp1);
		printf("Obtenido 1: %.20lf \n", obt1);
		printf("Esperado 2: %.20lf \n", esp2);
		printf("Obtenido 2: %.20lf \n", obt2);
		printf("Esperado 3: %.20lf \n", esp3);
		printf("Obtenido 3: %.20lf \n", obt3);
		printf("Esperado 4: %.20lf \n", esp4);
		printf("Obtenido 4: %.20lf \n", obt4);
		printf("Esperado 5: %.20lf \n", esp5);
		printf("Obtenido 5: %.20lf \n", obt5);
	}

	assert(fabs(esp1 - obt1) < 10e-12);
	assert(fabs(esp2 - obt2) < 10e-12);
	assert(fabs(esp3 - obt3) < 10e-12);
	assert(fabs(esp4 - obt4) < 10e-12);
	assert(fabs(esp5 - obt5) < 10e-12);
	printf("timediff1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testTimediff2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función timediff
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testTimediff2(bool verbose)
{
	double n1 = -0.141248008109364;
	double n2 = 34.0;
	double esp1 = -34.1412480081094;
	double esp2 = -15.0;
	double esp3 = -15.1412480081094;
	double esp4 = 66.184;
	double esp5 = 15.0;
	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;

	timediff(n1, n2, &obt1, &obt2, &obt3, &obt4, &obt5);

	if(verbose)
	{
		printf("timediff2:\n");
		printf("Esperado 1: %.20lf \n", esp1);
		printf("Obtenido 1: %.20lf \n", obt1);
		printf("Esperado 2: %.20lf \n", esp2);
		printf("Obtenido 2: %.20lf \n", obt2);
		printf("Esperado 3: %.20lf \n", esp3);
		printf("Obtenido 3: %.20lf \n", obt3);
		printf("Esperado 4: %.20lf \n", esp4);
		printf("Obtenido 4: %.20lf \n", obt4);
		printf("Esperado 5: %.20lf \n", esp5);
		printf("Obtenido 5: %.20lf \n", obt5);
	}

	assert(fabs(esp1 - obt1) < 10e-12);
	assert(fabs(esp2 - obt2) < 10e-12);
	assert(fabs(esp3 - obt3) < 10e-12);
	assert(fabs(esp4 - obt4) < 10e-12);
	assert(fabs(esp5 - obt5) < 10e-12);
	printf("timediff2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testTimediff3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función timediff
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testTimediff3(bool verbose)
{
	double n1 = 0.161894409590604;
	double n2 = 33.0;
	double esp1 = -32.8381055904094;
	double esp2 = -14.0;
	double esp3 = -13.8381055904094;
	double esp4 = 65.184;
	double esp5 = 14.0;
	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;

	timediff(n1, n2, &obt1, &obt2, &obt3, &obt4, &obt5);

	if(verbose)
	{
		printf("timediff3:\n");
		printf("Esperado 1: %.20lf \n", esp1);
		printf("Obtenido 1: %.20lf \n", obt1);
		printf("Esperado 2: %.20lf \n", esp2);
		printf("Obtenido 2: %.20lf \n", obt2);
		printf("Esperado 3: %.20lf \n", esp3);
		printf("Obtenido 3: %.20lf \n", obt3);
		printf("Esperado 4: %.20lf \n", esp4);
		printf("Obtenido 4: %.20lf \n", obt4);
		printf("Esperado 5: %.20lf \n", esp5);
		printf("Obtenido 5: %.20lf \n", obt5);
	}

	assert(fabs(esp1 - obt1) < 10e-12);
	assert(fabs(esp2 - obt2) < 10e-12);
	assert(fabs(esp3 - obt3) < 10e-12);
	assert(fabs(esp4 - obt4) < 10e-12);
	assert(fabs(esp5 - obt5) < 10e-12);
	printf("timediff3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testIERS1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función IERS
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testIERS1(bool verbose)
{
	double *eopdata[13];
	leerFichero(eopdata);
	int n1 = 20026;
	double n2 = 54977.6669036457;
	char c = 'l';
	double esp1 = 0.258022690875596;
	double esp2 = 34.0;
	double esp3 = 7.5789200806793e-08;
	double esp4 = 2.56777193042581e-06;
	double esp5 = -2.91335538497448e-07;
	double esp6 = -4.54223178651006e-08;
	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;
	double obt6;

	IERS(eopdata, n1, n2, c, &obt1, &obt2, &obt3, &obt4, &obt5, &obt6);

	if(verbose)
	{
		printf("IERS1:\n");
		printf("Esperado 1: %.20lf \n", esp1);
		printf("Obtenido 1: %.20lf \n", obt1);
		printf("Esperado 2: %.20lf \n", esp2);
		printf("Obtenido 2: %.20lf \n", obt2);
		printf("Esperado 3: %.20lf \n", esp3);
		printf("Obtenido 3: %.20lf \n", obt3);
		printf("Esperado 4: %.20lf \n", esp4);
		printf("Obtenido 4: %.20lf \n", obt4);
		printf("Esperado 5: %.20lf \n", esp5);
		printf("Obtenido 5: %.20lf \n", obt5);
		printf("Esperado 6: %.20lf \n", esp6);
		printf("Obtenido 6: %.20lf \n", obt6);
	}

	assert(fabs(esp1 - obt1) < 10e-12);
	assert(fabs(esp2 - obt2) < 10e-12);
	assert(fabs(esp3 - obt3) < 10e-12);
	assert(fabs(esp4 - obt4) < 10e-12);
	assert(fabs(esp5 - obt5) < 10e-12);
	assert(fabs(esp6 - obt6) < 10e-12);
	printf("IERS1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testIERS2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función IERS
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testIERS2(bool verbose)
{
	double *eopdata[13];
	leerFichero(eopdata);
	int n1 = 20026;
	double n2 = 53989.198426771;
	char c = 'l';
	double esp1 = 0.161894409590604;
	double esp2 = 33.0;
	double esp3 = 3.41489253490265e-07;
	double esp4 = 1.23201492955798e-06;
	double esp5 = -3.25133150333998e-07;
	double esp6 = -2.74651440721866e-08;
	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;
	double obt6;

	IERS(eopdata, n1, n2, c, &obt1, &obt2, &obt3, &obt4, &obt5, &obt6);

	if(verbose)
	{
		printf("IERS2:\n");
		printf("Esperado 1: %.20lf \n", esp1);
		printf("Obtenido 1: %.20lf \n", obt1);
		printf("Esperado 2: %.20lf \n", esp2);
		printf("Obtenido 2: %.20lf \n", obt2);
		printf("Esperado 3: %.20lf \n", esp3);
		printf("Obtenido 3: %.20lf \n", obt3);
		printf("Esperado 4: %.20lf \n", esp4);
		printf("Obtenido 4: %.20lf \n", obt4);
		printf("Esperado 5: %.20lf \n", esp5);
		printf("Obtenido 5: %.20lf \n", obt5);
		printf("Esperado 6: %.20lf \n", esp6);
		printf("Obtenido 6: %.20lf \n", obt6);
	}

	assert(fabs(esp1 - obt1) < 10e-12);
	assert(fabs(esp2 - obt2) < 10e-12);
	assert(fabs(esp3 - obt3) < 10e-12);
	assert(fabs(esp4 - obt4) < 10e-12);
	assert(fabs(esp5 - obt5) < 10e-12);
	assert(fabs(esp6 - obt6) < 10e-12);
	printf("IERS2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testIERS3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función IERS
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testIERS3(bool verbose)
{
	double *eopdata[13];
	leerFichero(eopdata);
	int n1 = 20026;
	double n2 = 55565.9044073611;
	char c = 'l';
	double esp1 = -0.141329177689613;
	double esp2 = 34.0;
	double esp3 = 5.72229229392844e-07;
	double esp4 = 9.71655423937364e-07;
	double esp5 = -3.22818789921021e-07;
	double esp6 = -3.03101990895201e-08;
	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;
	double obt6;

	IERS(eopdata, n1, n2, c, &obt1, &obt2, &obt3, &obt4, &obt5, &obt6);

	if(verbose)
	{
		printf("IERS3:\n");
		printf("Esperado 1: %.20lf \n", esp1);
		printf("Obtenido 1: %.20lf \n", obt1);
		printf("Esperado 2: %.20lf \n", esp2);
		printf("Obtenido 2: %.20lf \n", obt2);
		printf("Esperado 3: %.20lf \n", esp3);
		printf("Obtenido 3: %.20lf \n", obt3);
		printf("Esperado 4: %.20lf \n", esp4);
		printf("Obtenido 4: %.20lf \n", obt4);
		printf("Esperado 5: %.20lf \n", esp5);
		printf("Obtenido 5: %.20lf \n", obt5);
		printf("Esperado 6: %.20lf \n", esp6);
		printf("Obtenido 6: %.20lf \n", obt6);
	}

	assert(fabs(esp1 - obt1) < 10e-12);
	assert(fabs(esp2 - obt2) < 10e-12);
	assert(fabs(esp3 - obt3) < 10e-12);
	assert(fabs(esp4 - obt4) < 10e-12);
	assert(fabs(esp5 - obt5) < 10e-12);
	assert(fabs(esp6 - obt6) < 10e-12);
	printf("IERS3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGstime1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función gstime
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGstime1(bool verbose)
{
	double n = 54977.6669066321;
	double esperado = 2.17192854407046;
	double obtenido = gstime(n);

	if(verbose)
	{
		printf("gstime1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-8);
	printf("gstime1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGstime2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función gstime
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGstime2(bool verbose)
{
	double n = 55565.9044057253;
	double esperado = 1.2171566871639;
	double obtenido = gstime(n);

	if(verbose)
	{
		printf("gstime2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-8);
	printf("gstime2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGstime3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función gstime
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGstime3(bool verbose)
{
	double n = 54332.4861092138;
	double esperado = 2.50337730618796;
	double obtenido = gstime(n);

	if(verbose)
	{
		printf("gstime3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-8);
	printf("gstime3 superado!\n");
}

//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
 * Función principal, realiza todas las comprobaciones.
 */
//------------------------------------------------------------------------------
int main(){

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

	// Test gstime
	printf("Probando gstime!\n");
	testGstime1(false);
	testGstime2(false);
	testGstime3(false);
	printf("gstime finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
