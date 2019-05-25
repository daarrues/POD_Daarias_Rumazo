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

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
