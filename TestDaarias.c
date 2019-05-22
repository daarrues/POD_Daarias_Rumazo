#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "lambert_gooding.h"

typedef int bool;
#define true 1
#define false 0

//------------------------------------------------------------------------------
//  void testD8rt1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función d8rt
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testD8rt1(bool verbose)
{
	double n = 0.1373209352524762;
	double esperado = 0.7802200275631945;
	double obtenido = d8rt(n);

	if(verbose)
	{
		printf("d8rt1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("d8rt1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testD8rt2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función d8rt
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testD8rt2(bool verbose)
{
	double n = 0.005544910621842289;
	double esperado = 0.5223803340446431;
	double obtenido = d8rt(n);

	if(verbose)
	{
		printf("d8rt2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("d8rt2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testD8rt3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función d8rt
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testD8rt3(bool verbose)
{
	double n = 0.006959569773385268;
	double esperado = 0.537431105725774;
	double obtenido = d8rt(n);

	if(verbose)
	{
		printf("d8rt3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("d8rt3 superado!\n");
}

//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
 * Función principal, realiza todas las comprobaciones.
 */
//------------------------------------------------------------------------------
int main(){

	// Test d8rt
	printf("Probando d8rt!\n");
	testD8rt1(true);
	testD8rt2(true);
	testD8rt3(true);
	printf("d8rt finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
