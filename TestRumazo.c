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
#include "GHAMatrix.h"
#include "PoleMatrix.h"
#include "NutMatrix.h"
#include "PrecMatrix.h"
#include "doubler.h"
#include "rv2coe.h"

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
	double n = 54977.66766966425300000000;
	double esperado = 0.40907147055862825000;
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
	double n = 55565.90517337957900000000;
	double esperado = 0.40906781525909947000;
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
	double n = 54332.48686555545300000000;
	double esperado = 0.40907547970278441000;
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
	double n = 54977.66766966425300000000;
	double esperado1 = 0.00006486933873938495;
	double esperado2 = 0.00002230513367072657;
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
	double n = 55565.90517337957900000000;
	double esperado1 = 0.00008742355142398978;
	double esperado2 = -0.00000075329307704148;
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
	double n = 54332.48686555545300000000;
	double esperado1 = 0.00003520723806122855;
	double esperado2 = 0.00004182601666776308;
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
	double n = 54977.66766966426100000000;
	double esperado = 0.00005951700514220537;
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
	double n = 55565.90517337957900000000;
	double esperado = 0.00008021040920233633;
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
	double n = 54332.48686555545300000000;
	double esperado = 0.00003230225199490219;
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
	double n = 10.34566369140563100000;
	double esperado = 0.34566369140563147000;
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
	double n = 12.19370373739637800000;
	double esperado = 0.19370373739637792000;
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
	double n = 8.39841973164080890000;
	double esperado = 0.39841973164080891000;
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
	double n = 54977.66690663211200000000;
	double esperado = 2.17186902706532250000;
	double obtenido = gmst(n);

	if(verbose)
	{
		printf("gmst1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
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
	double n = 55565.90440572530500000000;
	double esperado = 1.21707647675469470000;
	double obtenido = gmst(n);

	if(verbose)
	{
		printf("gmst2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
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
	double n = 54332.48610921379400000000;
	double esperado = 2.50334500393596440000;
	double obtenido = gmst(n);

	if(verbose)
	{
		printf("gmst3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
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
	double n1 = 0.25802269087559637000;
	double n2 = 34.00000000000000000000;
	double esp1 = -33.74197730912440100000;
	double esp2 = -15.00000000000000000000;
	double esp3 = -14.74197730912440100000;
	double esp4 = 66.18399999999999700000;
	double esp5 = 15.00000000000000000000;
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
	double n1 = -0.14132917768961278000;
	double n2 = 34.00000000000000000000;
	double esp1 = -34.14132917768961300000;
	double esp2 = -15.00000000000000000000;
	double esp3 = -15.14132917768961300000;
	double esp4 = 66.18399999999999700000;
	double esp5 = 15.00000000000000000000;
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
	double n1 = -0.16391912638890360000;
	double n2 = 33.00000000000000000000;
	double esp1 = -33.16391912638890500000;
	double esp2 = -14.00000000000000000000;
	double esp3 = -14.16391912638890500000;
	double esp4 = 65.18399999999999700000;
	double esp5 = 14.00000000000000000000;
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
	double n2 = 54977.66690364573200000000;
	char c = 'l';
	double esp1 = 0.25802269087559637000;
	double esp2 = 34.00000000000000000000;
	double esp3 = 0.00000007578920080679;
	double esp4 = 0.00000256777193042581;
	double esp5 = -0.00000029133553849745;
	double esp6 = -0.00000004542231786510;
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
	double n2 = 55565.90440736105700000000;
	char c = 'l';
	double esp1 = -0.14132917768961278000;
	double esp2 = 34.00000000000000000000;
	double esp3 = 0.00000057222922939284;
	double esp4 = 0.00000097165542393736;
	double esp5 = -0.00000032281878992102;
	double esp6 = -0.00000003031019908952;
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
	double n2 = 54332.48611111100800000000;
	char c = 'l';
	double esp1 = -0.16391912638890360000;
	double esp2 = 33.00000000000000000000;
	double esp3 = 0.00000101984436330847;
	double esp4 = 0.00000138460955806667;
	double esp5 = -0.00000032922215709019;
	double esp6 = -0.00000003203992214460;
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
//  void testGast1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función gast
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGast1(bool verbose)
{
	double n = 54977.66690663211200000000;
	double esperado = 2.17192854407046450000;
	double obtenido = gast(n);

	if(verbose)
	{
		printf("gast1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("gast1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGast2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función gast
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGast2(bool verbose)
{
	double n = 55565.90440572530500000000;
	double esperado = 1.21715668716389700000;
	double obtenido = gast(n);

	if(verbose)
	{
		printf("gast2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("gast2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGast3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función gast
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGast3(bool verbose)
{
	double n = 54332.48610921379400000000;
	double esperado = 2.50337730618795940000;
	double obtenido = gast(n);

	if(verbose)
	{
		printf("gast3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("gast3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGHAMatrix1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función GHAMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGHAMatrix1(bool verbose)
{
	double date = 54977.66690663211200000000;
	double esperado[3][3] = {
		{
			-0.56557657052466759000,
			0.82469578807797705000,
			0.00000000000000000000
		},
		{
			-0.82469578807797705000,
			-0.56557657052466759000,
			0.00000000000000000000
		},
		{
			0.00000000000000000000,
			0.00000000000000000000,
			1.00000000000000000000
		}
	};

	double obtenido[3][3];
	GHAMatrix(date, obtenido);

	if(verbose)
	{
		printf("GHAMatrix1:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("GHAMatrix1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGHAMatrix2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función GHAMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGHAMatrix2(bool verbose)
{
	double date = 55565.90440572530500000000;
	double esperado[3][3] = {
		{
			0.34631450688377857000,
			0.93811846923607967000,
			0.00000000000000000000
		},
		{
			-0.93811846923607967000,
			0.34631450688377857000,
			0.00000000000000000000
		},
		{
			0.00000000000000000000,
			0.00000000000000000000,
			1.00000000000000000000
		}
	};

	double obtenido[3][3];
	GHAMatrix(date, obtenido);

	if(verbose)
	{
		printf("GHAMatrix2:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("GHAMatrix2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGHAMatrix3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función GHAMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGHAMatrix3(bool verbose)
{
	double date = 54332.48610921379400000000;
	double esperado[3][3] = {
		{
			-0.80316026638348281000,
			0.59576302881499199000,
			0.00000000000000000000
		},
		{
			-0.59576302881499199000,
			-0.80316026638348281000,
			0.00000000000000000000
		},
		{
			0.00000000000000000000,
			0.00000000000000000000,
			1.00000000000000000000
		}
	};

	double obtenido[3][3];
	GHAMatrix(date, obtenido);

	if(verbose)
	{
		printf("GHAMatrix3:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("GHAMatrix3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPoleMatrix1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función PoleMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPoleMatrix1(bool verbose)
{
	double xp = 0.00000007578920080679;
	double yp = 0.00000256777193042581;
	double esperado[3][3] = {
		{
			0.99999999999999711000,
			0.00000000000019460938,
			0.00000007578920080654
		},
		{
			0.00000000000000000000,
			0.99999999999670330000,
			-0.00000256777193042298
		},
		{
			-0.00000007578920080679,
			0.00000256777193042298,
			0.99999999999670042000
		}
	};

	double obtenido[3][3];
	PoleMatrix(xp, yp, obtenido);

	if(verbose)
	{
		printf("PoleMatrix1:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("PoleMatrix1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPoleMatrix2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función PoleMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPoleMatrix2(bool verbose)
{
	double xp = 0.00000057222922939284;
	double yp = 0.00000097165542393736;
	double esperado[3][3] = {
		{
			0.99999999999983624000,
			0.00000000000055600963,
			0.00000057222922939254
		},
		{
			0.00000000000000000000,
			0.99999999999952793000,
			-0.00000097165542393721
		},
		{
			-0.00000057222922939281,
			0.00000097165542393705,
			0.99999999999936418000
		}
	};

	double obtenido[3][3];
	PoleMatrix(xp, yp, obtenido);

	if(verbose)
	{
		printf("PoleMatrix2:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("PoleMatrix2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPoleMatrix3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función PoleMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPoleMatrix3(bool verbose)
{
	double xp = 0.00000101984436330847;
	double yp = 0.00000138460955806667;
	double esperado[3][3] = {
		{
			0.99999999999947997000,
			0.00000000000141208625,
			0.00000101984436330732
		},
		{
			0.00000000000000000000,
			0.99999999999904143000,
			-0.00000138460955806623
		},
		{
			-0.00000101984436330829,
			0.00000138460955806551,
			0.99999999999852140000
		}
	};

	double obtenido[3][3];
	PoleMatrix(xp, yp, obtenido);

	if(verbose)
	{
		printf("PoleMatrix3:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("PoleMatrix3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutMatrix1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función NutMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNutMatrix1(bool verbose)
{
	double date = 54977.66766966425300000000;
	double esperado[3][3] = {
		{
			0.99999999789598448000,
			-0.00005951700510045734,
			-0.00002580227134293989
		},
		{
			0.00005951642956254065,
			0.99999999798012063000,
			-0.00002230590149845583
		},
		{
			0.00002580359887127566,
			0.00002230436579244361,
			0.99999999941834472000
		}
	};

	double obtenido[3][3];
	NutMatrix(date, obtenido);

	if(verbose)
	{
		printf("NutMatrix1:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("NutMatrix1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutMatrix2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función NutMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNutMatrix2(bool verbose)
{
	double date = 55565.90517337957900000000;
	double esperado[3][3] = {
		{
			0.99999999617856128000,
			-0.00008021040910016327,
			-0.00003477308723849868
		},
		{
			0.00008021043529446639,
			0.99999999678286022000,
			0.00000075189849470954
		},
		{
			0.00003477302681654293,
			-0.00000075468765634400,
			0.99999999939513351000
		}
	};

	double obtenido[3][3];
	NutMatrix(date, obtenido);

	if(verbose)
	{
		printf("NutMatrix2:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("NutMatrix2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutMatrix3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función NutMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNutMatrix3(bool verbose)
{
	double date = 54332.48686555545300000000;
	double esperado[3][3] = {
		{
			0.99999999938022521000,
			-0.00003230225198822882,
			-0.00001400407540242579
		},
		{
			0.00003230166622528275,
			0.99999999860358380000,
			-0.00004182624283308504
		},
		{
			0.00001400542646470598,
			0.00004182579045223188,
			0.99999999902722569000
		}
	};

	double obtenido[3][3];
	NutMatrix(date, obtenido);

	if(verbose)
	{
		printf("NutMatrix3:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("NutMatrix3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPrecMatrix1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función PrecMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPrecMatrix1(bool verbose)
{
	double date1 = 51544.50000000000000000000;
	double date2 = 54977.66766966425300000000;
	double esperado[3][3] = {
		{
			0.99999737380232945000,
			-0.00210194819368334690,
			-0.00091334672251535971
		},
		{
			0.00210194819366918240,
			0.99999779090399488000,
			-0.00000095992051556750
		},
		{
			0.00091334672254795744,
			-0.00000095988949895835,
			0.99999958289833457000
		}
	};

	double obtenido[3][3];
	PrecMatrix(date1, date2, obtenido);

	if(verbose)
	{
		printf("PrecMatrix1:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("PrecMatrix1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPrecMatrix2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función PrecMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPrecMatrix2(bool verbose)
{
	double date1 = 51544.50000000000000000000;
	double date2 = 55565.90517337957900000000;
	double esperado[3][3] = {
		{
			0.99999639673612872000,
			-0.00246210631618178730,
			-0.00106983514929104710
		},
		{
			0.00246210631615512330,
			0.99999696901078317000,
			-0.00000131705123572390
		},
		{
			0.00106983514935241200,
			-0.00000131700138827348,
			0.99999942772534556000
		}
	};

	double obtenido[3][3];
	PrecMatrix(date1, date2, obtenido);

	if(verbose)
	{
		printf("PrecMatrix2:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("PrecMatrix2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPrecMatrix3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función PrecMatrix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPrecMatrix3(bool verbose)
{
	double date1 = 51544.50000000000000000000;
	double date2 = 54332.48686555545300000000;
	double esperado[3][3] = {
		{
			0.99999826812913417000,
			-0.00170692926348680160,
			-0.00074170831295002732
		},
		{
			0.00170692926348064160,
			0.99999854319498327000,
			-0.00000063303066555599
		},
		{
			0.00074170831296420381,
			-0.00000063301405511264,
			0.99999972493415101000
		}
	};

	double obtenido[3][3];
	PrecMatrix(date1, date2, obtenido);

	if(verbose)
	{
		printf("PrecMatrix3:\n");
	}
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(esperado[i][j] - obtenido[i][j]) < 10e-12);
		}
	}

	printf("PrecMatrix3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testDoubler1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función doubler
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testDoubler1(bool verbose)
{
	double n1 = 5972180.93003293410000000000;
	double n2 = 6395944.28126916850000000000;
	double n3 = 6369760.82916144930000000000;
	double n4 = 6369760.82916145030000000000;
	double n5 = 12820055.36999999900000000000;
	double n6 = 13457869.06999999800000000000;
	double v1[3] = {
		0.63388609516515371000,
		-0.77334037977825631000,
		0.01153582943251443100
	};
	double v2[3] = {
		0.93553982556964899000,
		-0.01648308182240541800,
		-0.35283642497161005000
	};
	double v3[3] = {
		0.59638436830343833000,
		0.53607267396708524000,
		-0.59745441120564824000
	};
	double v4[3] = {
		4950990.33826459660000000000,
		256563.11626038095000000000,
		3999465.34658133240000000000
	};
	double v5[3] = {
		4935037.85913035740000000000,
		472703.32020261511000000000,
		3999475.70573182080000000000
	};
	double v6[3] = {
		4909646.95198535550000000000,
		687938.93691575713000000000,
		3999494.94894738870000000000
	};
	double n7 = -600.00000447034836000000;
	double n8 = 600.00000447034836000000;
	char c = 'y';

	double esp1[3] = {
		13430456.18911721200000000000,
		323024.30962337909000000000,
		795450.66022862401000000000
	};
	double esp2[3] = {
		11608410.68523451300000000000,
		6709264.07153013260000000000,
		-2711287.81928542350000000000
	};
	double esp3 = 933.40653066845334000000;
	double esp4 = -879.39690150175784000000;
	double esp5 = 1282.41438773331400000000;
	double esp6 = 12820055.36999999900000000000;
	double esp7 = 13457869.07000000000000000000;
	double esp8 = 12459660.56374719900000000000;
	double esp9 = 0.61504410921129493000;

	double obt1[3];
	double obt2[3];
	double obt3;
	double obt4;
	double obt5;
	double obt6;
	double obt7;
	double obt8;
	double obt9;

	doubler(n1, n2, n3, n4, n5, n6, v1, v2, v3, v4, v5, v6, n7, n8, c,
	        obt1, obt2, &obt3, &obt4, &obt5, &obt6, &obt7, &obt8, &obt9);

	if(verbose)
	{
		printf("doubler1:\n");
		for(int i = 0; i < 3; i++)
		{
			printf("Esperado 1: %.20lf\n", esp1[i]);
			printf("Obtenido 1: %.20lf\n", obt1[i]);
			printf("Esperado 2: %.20lf\n", esp2[i]);
			printf("Obtenido 2: %.20lf\n", obt2[i]);
		}
		printf("Esperado 3: %.20lf\n", esp3);
		printf("Obtenido 3: %.20lf\n", obt3);
		printf("Esperado 4: %.20lf\n", esp4);
		printf("Obtenido 4: %.20lf\n", obt4);
		printf("Esperado 5: %.20lf\n", esp5);
		printf("Obtenido 5: %.20lf\n", obt5);
		printf("Esperado 6: %.20lf\n", esp6);
		printf("Obtenido 6: %.20lf\n", obt6);
		printf("Esperado 7: %.20lf\n", esp7);
		printf("Obtenido 7: %.20lf\n", obt7);
		printf("Esperado 8: %.20lf\n", esp8);
		printf("Obtenido 8: %.20lf\n", obt8);
		printf("Esperado 9: %.20lf\n", esp9);
		printf("Obtenido 9: %.20lf\n", obt9);
	}

	for(int i = 0; i < 3; i++)
	{
		assert(fabs(esp1[i] - obt1[i]) < 10e-6);
		assert(fabs(esp2[i] - obt2[i]) < 10e-6);
	}
	assert(fabs(esp3 - obt3) < 10e-6);
	assert(fabs(esp4 - obt4) < 10e-6);
	assert(fabs(esp5 - obt5) < 10e-6);
	assert(fabs(esp6 - obt6) < 10e-6);
	assert(fabs(esp7 - obt7) < 10e-6);
	assert(fabs(esp8 - obt8) < 10e-6);
	assert(fabs(esp9 - obt9) < 10e-6);

	printf("doubler1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testDoubler2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función doubler
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testDoubler2(bool verbose)
{
	double n1 = 7515275.95733092820000000000;
	double n2 = 7515565.47400444190000000000;
	double n3 = 6372639.11744251940000000000;
	double n4 = 6372639.11744252030000000000;
	double n5 = 12820055.36999999900000000000;
	double n6 = 13457869.06999999800000000000;
	double v1[3] = {
		0.94986053669058212000,
		0.29900279622503534000,
		-0.09144555039743179600
	};
	double v2[3] = {
		0.94309959610195582000,
		0.31969674480435684000,
		-0.09141741187446586500
	};
	double v3[3] = {
		0.93588812100465057000,
		0.34023783543213759000,
		-0.09138730930954264100
	};
	double v4[3] = {
		4991697.66717447340000000000,
		-2304791.72214133110000000000,
		3222020.92454910930000000000
	};
	double v5[3] = {
		5040923.41475751530000000000,
		-2195094.88576021140000000000,
		3221983.68869511690000000000
	};
	double v6[3] = {
		5087737.81584266200000000000,
		-2084347.48894624690000000000,
		3221947.97967418980000000000
	};
	double n7 = -299.99998211860657000000;
	double n8 = 300.00002235174179000000;
	char c = 'y';

	double esp1[3] = {
		13224229.70471591500000000000,
		578924.30950820725000000000,
		2428751.79177130620000000000
	};
	double esp2[3] = {
		13948424.03334139100000000000,
		1136914.36001161720000000000,
		2356722.47474513020000000000
	};
	double esp3 = -254.05855596709932000000;
	double esp4 = 248.34817672512381000000;
	double esp5 = 355.27815404663016000000;
	double esp6 = 12820055.37000000100000000000;
	double esp7 = 13457869.07000000000000000000;
	double esp8 = -1547696.59548539830000000000;
	double esp9 = 0.05997878907998426000;

	double obt1[3];
	double obt2[3];
	double obt3;
	double obt4;
	double obt5;
	double obt6;
	double obt7;
	double obt8;
	double obt9;

	doubler(n1, n2, n3, n4, n5, n6, v1, v2, v3, v4, v5, v6, n7, n8, c,
	        obt1, obt2, &obt3, &obt4, &obt5, &obt6, &obt7, &obt8, &obt9);

	if(verbose)
	{
		printf("doubler2:\n");
		for(int i = 0; i < 3; i++)
		{
			printf("Esperado 1: %.20lf\n", esp1[i]);
			printf("Obtenido 1: %.20lf\n", obt1[i]);
			printf("Esperado 2: %.20lf\n", esp2[i]);
			printf("Obtenido 2: %.20lf\n", obt2[i]);
		}
		printf("Esperado 3: %.20lf\n", esp3);
		printf("Obtenido 3: %.20lf\n", obt3);
		printf("Esperado 4: %.20lf\n", esp4);
		printf("Obtenido 4: %.20lf\n", obt4);
		printf("Esperado 5: %.20lf\n", esp5);
		printf("Obtenido 5: %.20lf\n", obt5);
		printf("Esperado 6: %.20lf\n", esp6);
		printf("Obtenido 6: %.20lf\n", obt6);
		printf("Esperado 7: %.20lf\n", esp7);
		printf("Obtenido 7: %.20lf\n", obt7);
		printf("Esperado 8: %.20lf\n", esp8);
		printf("Obtenido 8: %.20lf\n", obt8);
		printf("Esperado 9: %.20lf\n", esp9);
		printf("Obtenido 9: %.20lf\n", obt9);
	}

	for(int i = 0; i < 3; i++)
	{
		assert(fabs(esp1[i] - obt1[i]) < 10e-6);
		assert(fabs(esp2[i] - obt2[i]) < 10e-6);
	}
	assert(fabs(esp3 - obt3) < 10e-6);
	assert(fabs(esp4 - obt4) < 10e-6);
	assert(fabs(esp5 - obt5) < 10e-6);
	assert(fabs(esp6 - obt6) < 10e-6);
	assert(fabs(esp7 - obt7) < 10e-6);
	assert(fabs(esp8 - obt8) < 10e-6);
	assert(fabs(esp9 - obt9) < 10e-6);

	printf("doubler2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testDoubler3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función doubler
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testDoubler3(bool verbose)
{
	double n1 = 10215883.52101078600000000000;
	double n2 = 12282028.53736285900000000000;
	double n3 = 6371344.85231744870000000000;
	double n4 = 6371344.85231744960000000000;
	double n5 = 12820055.36999999900000000000;
	double n6 = 13457869.06999999800000000000;
	double v1[3] = {
		0.95388721001727206000,
		-0.00694753042371069540,
		0.30008485864247841000
	};
	double v2[3] = {
		0.45956638535656152000,
		0.65860637495401231000,
		0.59584929329507408000
	};
	double v3[3] = {
		-0.68285277220415785000,
		0.70016805580162622000,
		0.20851087532324863000
	};
	double v4[3] = {
		4092160.77507185750000000000,
		2689587.44262100380000000000,
		4076073.45451609280000000000
	};
	double v5[3] = {
		3970615.77626515740000000000,
		2865864.80883913630000000000,
		4076158.05408520900000000000
	};
	double v6[3] = {
		3561420.55410003290000000000,
		3360490.91297653410000000000,
		4076446.94406636290000000000
	};
	double n7 = -600.00000447034836000000;
	double n8 = 1800.00001341104510000000;
	char c = 'y';

	double esp1[3] = {
		7283791.29584485110000000000,
		7613989.53698643110000000000,
		8371844.93373465070000000000
	};
	double esp2[3] = {
		-2886203.32121872530000000000,
		9971608.94175766410000000000,
		6045245.71405340920000000000
	};
	double esp3 = 794.18467835950469000000;
	double esp4 = -509.21835069695499000000;
	double esp5 = 943.41540798712674000000;
	double esp6 = 12820055.36999999700000000000;
	double esp7 = 13457869.06999999800000000000;
	double esp8 = 10603469.17439548700000000000;
	double esp9 = 1.09115287002394970000;

	double obt1[3];
	double obt2[3];
	double obt3;
	double obt4;
	double obt5;
	double obt6;
	double obt7;
	double obt8;
	double obt9;

	doubler(n1, n2, n3, n4, n5, n6, v1, v2, v3, v4, v5, v6, n7, n8, c,
	        obt1, obt2, &obt3, &obt4, &obt5, &obt6, &obt7, &obt8, &obt9);

	if(verbose)
	{
		printf("doubler3:\n");
		for(int i = 0; i < 3; i++)
		{
			printf("Esperado 1: %.20lf\n", esp1[i]);
			printf("Obtenido 1: %.20lf\n", obt1[i]);
			printf("Esperado 2: %.20lf\n", esp2[i]);
			printf("Obtenido 2: %.20lf\n", obt2[i]);
		}
		printf("Esperado 3: %.20lf\n", esp3);
		printf("Obtenido 3: %.20lf\n", obt3);
		printf("Esperado 4: %.20lf\n", esp4);
		printf("Obtenido 4: %.20lf\n", obt4);
		printf("Esperado 5: %.20lf\n", esp5);
		printf("Obtenido 5: %.20lf\n", obt5);
		printf("Esperado 6: %.20lf\n", esp6);
		printf("Obtenido 6: %.20lf\n", obt6);
		printf("Esperado 7: %.20lf\n", esp7);
		printf("Obtenido 7: %.20lf\n", obt7);
		printf("Esperado 8: %.20lf\n", esp8);
		printf("Obtenido 8: %.20lf\n", obt8);
		printf("Esperado 9: %.20lf\n", esp9);
		printf("Obtenido 9: %.20lf\n", obt9);
	}

	for(int i = 0; i < 3; i++)
	{
		assert(fabs(esp1[i] - obt1[i]) < 10e-6);
		assert(fabs(esp2[i] - obt2[i]) < 10e-6);
	}
	assert(fabs(esp3 - obt3) < 10e-6);
	assert(fabs(esp4 - obt4) < 10e-6);
	assert(fabs(esp5 - obt5) < 10e-6);
	assert(fabs(esp6 - obt6) < 10e-6);
	assert(fabs(esp7 - obt7) < 10e-6);
	assert(fabs(esp8 - obt8) < 10e-6);
	assert(fabs(esp9 - obt9) < 10e-6);

	printf("doubler3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testRv2coe1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función rv2coe
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testRv2coe1(bool verbose)
{
	double v1[3] = {
		20435422.35215374500000000000,
		1070699.44671813420000000000,
		1012905.49143374850000000000
	};
	double v2[3] = {
		17.19646978075456200000,
		-2657.51027611591140000000,
		3738.38685080886900000000
	};

	double esp1  = 22151801.03597825400000000000;
	double esp2  = 22303727.64121050000000000000;
	double esp3  = 0.08253310617383882700;
	double esp4  = 2.18724254265762050000;
	double esp5  = 0.08742599164324867100;
	double esp6  = 6.16261219030337860000;
	double esp7  = 0.18120031105310727000;
	double esp8  = 0.15317433135474373000;
	double esp9  = 999999.09999999998000000000;
	double esp10 = 999999.09999999998000000000;
	double esp11 = 999999.09999999998000000000;

	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;
	double obt6;
	double obt7;
	double obt8;
	double obt9;
	double obt10;
	double obt11;

	rv2coe(v1, v2, &obt1, &obt2, &obt3, &obt4, &obt5, &obt6, &obt7, &obt8,
		      &obt9, &obt10, &obt11);

	if(verbose)
	{
		printf("rv2coe1:\n");
		printf("Esperado 1:  %.20lf\n", esp1);
		printf("Obtenido 1:  %.20lf\n", obt1);
		printf("Esperado 2:  %.20lf\n", esp2);
		printf("Obtenido 2:  %.20lf\n", obt2);
		printf("Esperado 3:  %.20lf\n", esp3);
		printf("Obtenido 3:  %.20lf\n", obt3);
		printf("Esperado 4:  %.20lf\n", esp4);
		printf("Obtenido 4:  %.20lf\n", obt4);
		printf("Esperado 5:  %.20lf\n", esp5);
		printf("Obtenido 5:  %.20lf\n", obt5);
		printf("Esperado 6:  %.20lf\n", esp6);
		printf("Obtenido 6:  %.20lf\n", obt6);
		printf("Esperado 7:  %.20lf\n", esp7);
		printf("Obtenido 7:  %.20lf\n", obt7);
		printf("Esperado 8:  %.20lf\n", esp8);
		printf("Obtenido 8:  %.20lf\n", obt8);
		printf("Esperado 9:  %.20lf\n", esp9);
		printf("Obtenido 9:  %.20lf\n", obt9);
		printf("Esperado 10: %.20lf\n", esp10);
		printf("Obtenido 10: %.20lf\n", obt10);
		printf("Esperado 11: %.20lf\n", esp11);
		printf("Obtenido 11: %.20lf\n", obt11);
	}

	assert(fabs(esp1 - obt1) < 10e-6);
	assert(fabs(esp2 - obt2) < 10e-6);
	assert(fabs(esp3 - obt3) < 10e-6);
	assert(fabs(esp4 - obt4) < 10e-6);
	assert(fabs(esp5 - obt5) < 10e-6);
	assert(fabs(esp6 - obt6) < 10e-6);
	assert(fabs(esp7 - obt7) < 10e-6);
	assert(fabs(esp8 - obt8) < 10e-6);
	assert(fabs(esp9 - obt9) < 10e-6);
	assert(fabs(esp10 - obt10) < 10e-6);
	assert(fabs(esp11 - obt11) < 10e-6);

	printf("rv2coe1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testRv2coe2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función rv2coe
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testRv2coe2(bool verbose)
{
	double v1[3] = {
		20456329.59102045000000000000,
		1074191.36683772060000000000,
		1009857.02167864280000000000
	};
	double v2[3] = {
		17.65947568968910300000,
		-2661.67627260588730000000,
		3743.57355524540890000000
	};

	double esp1  = 22261413.24800409000000000000;
	double esp2  = 22429994.93127170600000000000;
	double esp3  = 0.08669431216174891500;
	double esp4  = 2.18733488568994480000;
	double esp5  = 0.08740807153658009000;
	double esp6  = 6.17226628945409580000;
	double esp7  = 0.17130585226709599000;
	double esp8  = 0.14354279174665671000;
	double esp9  = 999999.09999999998000000000;
	double esp10 = 999999.09999999998000000000;
	double esp11 = 999999.09999999998000000000;

	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;
	double obt6;
	double obt7;
	double obt8;
	double obt9;
	double obt10;
	double obt11;

	rv2coe(v1, v2, &obt1, &obt2, &obt3, &obt4, &obt5, &obt6, &obt7, &obt8,
		      &obt9, &obt10, &obt11);

	if(verbose)
	{
		printf("rv2coe2:\n");
		printf("Esperado 1:  %.20lf\n", esp1);
		printf("Obtenido 1:  %.20lf\n", obt1);
		printf("Esperado 2:  %.20lf\n", esp2);
		printf("Obtenido 2:  %.20lf\n", obt2);
		printf("Esperado 3:  %.20lf\n", esp3);
		printf("Obtenido 3:  %.20lf\n", obt3);
		printf("Esperado 4:  %.20lf\n", esp4);
		printf("Obtenido 4:  %.20lf\n", obt4);
		printf("Esperado 5:  %.20lf\n", esp5);
		printf("Obtenido 5:  %.20lf\n", obt5);
		printf("Esperado 6:  %.20lf\n", esp6);
		printf("Obtenido 6:  %.20lf\n", obt6);
		printf("Esperado 7:  %.20lf\n", esp7);
		printf("Obtenido 7:  %.20lf\n", obt7);
		printf("Esperado 8:  %.20lf\n", esp8);
		printf("Obtenido 8:  %.20lf\n", obt8);
		printf("Esperado 9:  %.20lf\n", esp9);
		printf("Obtenido 9:  %.20lf\n", obt9);
		printf("Esperado 10: %.20lf\n", esp10);
		printf("Obtenido 10: %.20lf\n", obt10);
		printf("Esperado 11: %.20lf\n", esp11);
		printf("Obtenido 11: %.20lf\n", obt11);
	}

	assert(fabs(esp1 - obt1) < 10e-6);
	assert(fabs(esp2 - obt2) < 10e-6);
	assert(fabs(esp3 - obt3) < 10e-6);
	assert(fabs(esp4 - obt4) < 10e-6);
	assert(fabs(esp5 - obt5) < 10e-6);
	assert(fabs(esp6 - obt6) < 10e-6);
	assert(fabs(esp7 - obt7) < 10e-6);
	assert(fabs(esp8 - obt8) < 10e-6);
	assert(fabs(esp9 - obt9) < 10e-6);
	assert(fabs(esp10 - obt10) < 10e-6);
	assert(fabs(esp11 - obt11) < 10e-6);

	printf("rv2coe2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testRv2coe3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función rv2coe
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testRv2coe3(bool verbose)
{
	double v1[3] = {
		20418280.37423648300000000000,
		1067836.39923680880000000000,
		1015404.95114512560000000000
	};
	double v2[3] = {
		16.87979502296278200000,
		-2654.08002932726280000000,
		3734.12004615458320000000
	};

	double esp1  = 22062031.73590861600000000000;
	double esp2  = 22201041.68683477100000000000;
	double esp3  = 0.07912910777856066400;
	double esp4  = 2.18716682303263490000;
	double esp5  = 0.08744068138692351300;
	double esp6  = 6.15374268501359100000;
	double esp7  = 0.19026723745652355000;
	double esp8  = 0.16199787055311890000;
	double esp9  = 999999.09999999998000000000;
	double esp10 = 999999.09999999998000000000;
	double esp11 = 999999.09999999998000000000;

	double obt1;
	double obt2;
	double obt3;
	double obt4;
	double obt5;
	double obt6;
	double obt7;
	double obt8;
	double obt9;
	double obt10;
	double obt11;

	rv2coe(v1, v2, &obt1, &obt2, &obt3, &obt4, &obt5, &obt6, &obt7, &obt8,
		      &obt9, &obt10, &obt11);

	if(verbose)
	{
		printf("rv2coe3:\n");
		printf("Esperado 1:  %.20lf\n", esp1);
		printf("Obtenido 1:  %.20lf\n", obt1);
		printf("Esperado 2:  %.20lf\n", esp2);
		printf("Obtenido 2:  %.20lf\n", obt2);
		printf("Esperado 3:  %.20lf\n", esp3);
		printf("Obtenido 3:  %.20lf\n", obt3);
		printf("Esperado 4:  %.20lf\n", esp4);
		printf("Obtenido 4:  %.20lf\n", obt4);
		printf("Esperado 5:  %.20lf\n", esp5);
		printf("Obtenido 5:  %.20lf\n", obt5);
		printf("Esperado 6:  %.20lf\n", esp6);
		printf("Obtenido 6:  %.20lf\n", obt6);
		printf("Esperado 7:  %.20lf\n", esp7);
		printf("Obtenido 7:  %.20lf\n", obt7);
		printf("Esperado 8:  %.20lf\n", esp8);
		printf("Obtenido 8:  %.20lf\n", obt8);
		printf("Esperado 9:  %.20lf\n", esp9);
		printf("Obtenido 9:  %.20lf\n", obt9);
		printf("Esperado 10: %.20lf\n", esp10);
		printf("Obtenido 10: %.20lf\n", obt10);
		printf("Esperado 11: %.20lf\n", esp11);
		printf("Obtenido 11: %.20lf\n", obt11);
	}

	assert(fabs(esp1 - obt1) < 10e-6);
	assert(fabs(esp2 - obt2) < 10e-6);
	assert(fabs(esp3 - obt3) < 10e-6);
	assert(fabs(esp4 - obt4) < 10e-6);
	assert(fabs(esp5 - obt5) < 10e-6);
	assert(fabs(esp6 - obt6) < 10e-6);
	assert(fabs(esp7 - obt7) < 10e-6);
	assert(fabs(esp8 - obt8) < 10e-6);
	assert(fabs(esp9 - obt9) < 10e-6);
	assert(fabs(esp10 - obt10) < 10e-6);
	assert(fabs(esp11 - obt11) < 10e-6);

	printf("rv2coe3 superado!\n");
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

	// Test gast
	printf("Probando gast!\n");
	testGast1(false);
	testGast2(false);
	testGast3(false);
	printf("gast finalizado!\n\n");

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

	// Test doubler
	printf("Probando doubler!\n");
	testDoubler1(false);
	testDoubler2(false);
	testDoubler3(false);
	printf("doubler finalizado!\n\n");

	// Test rv2coe
	printf("Probando rv2coe!\n");
	testRv2coe1(false);
	testRv2coe2(false);
	testRv2coe3(false);
	printf("rv2coe finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
