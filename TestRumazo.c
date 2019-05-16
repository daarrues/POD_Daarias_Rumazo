#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "MeanObliquity.h"
#include "NutAngles.h"

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
 * Comprobación 1 de la función testNutAngles
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
		printf("testNutAngles1:\n");
		printf("Esperado 1: %.20lf \n", esperado1);
		printf("Obtenido 1: %.20lf \n", obtenido1);
		printf("Esperado 2: %.20lf \n", esperado2);
		printf("Obtenido 2: %.20lf \n", obtenido2);
	}

	assert(fabs(esperado1 - obtenido1) < 10e-12);
	assert(fabs(esperado2 - obtenido2) < 10e-12);
	printf("testNutAngles1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutAngles2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función testNutAngles
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
		printf("testNutAngles2:\n");
		printf("Esperado 1: %.20lf \n", esperado1);
		printf("Obtenido 1: %.20lf \n", obtenido1);
		printf("Esperado 2: %.20lf \n", esperado2);
		printf("Obtenido 2: %.20lf \n", obtenido2);
	}

	assert(fabs(esperado1 - obtenido1) < 10e-12);
	assert(fabs(esperado2 - obtenido2) < 10e-12);
	printf("testNutAngles2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNutAngles3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función testNutAngles
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
		printf("testNutAngles3:\n");
		printf("Esperado 1: %.20lf \n", esperado1);
		printf("Obtenido 1: %.20lf \n", obtenido1);
		printf("Esperado 2: %.20lf \n", esperado2);
		printf("Obtenido 2: %.20lf \n", obtenido2);
	}

	assert(fabs(esperado1 - obtenido1) < 10e-12);
	assert(fabs(esperado2 - obtenido2) < 10e-12);
	printf("testNutAngles3 superado!\n");
}

//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
 * Función principal, realiza todas las comprobaciones.
 */
//------------------------------------------------------------------------------
int main(){

	// Test testMeanObliquity
	printf("Probando MeanObliquity!\n");
	testMeanObliquity1(false);
	testMeanObliquity2(false);
	testMeanObliquity3(false);
	printf("MeanObliquity finalizado!\n\n");

	// Test testNutAngles
	printf("Probando NutAngles!\n");
	testNutAngles1(false);
	testNutAngles2(false);
	testNutAngles3(false);
	printf("NutAngles finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
