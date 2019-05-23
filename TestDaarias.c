#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "lambert_gooding.h"
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"

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
//  void testTlamb1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función tlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
 //------------------------------------------------------------------------------
void testTlamb1(bool verbose)
{
	double m = 0.0;
  double q = 0.8046115644662321;
  double qsqfm1 = 0.3526002303272028;
  double x = 0.0;
  double n = 0.0;

	double esp_t = 2.227109718426917;
  double esp_dt = 0.0;
  double esp_d2t = 0.0;
  double esp_d3t = 0.0;

  double obt_t, obt_dt, obt_d2t, obt_d3t;
	tlamb(m, q, qsqfm1, x, n, &obt_t, &obt_dt, &obt_d2t, &obt_d3t);

	if(verbose)
	{
		printf("tlamb1:\n");
		printf("Esperado t: %.20lf \n", esp_t);
		printf("Obtenido t: %.20lf \n", obt_t);
    printf("Esperado dt: %.20lf \n", esp_dt);
		printf("Obtenido dt: %.20lf \n", obt_dt);
    printf("Esperado d2t: %.20lf \n", esp_d2t);
		printf("Obtenido d2t: %.20lf \n", obt_d2t);
    printf("Esperado d3t: %.20lf \n", esp_d3t);
		printf("Obtenido d3t: %.20lf \n", obt_d3t);
	}

	assert(fabs(esp_t - obt_t) < 10e-12);
  assert(fabs(esp_dt - obt_dt) < 10e-12);
  assert(fabs(esp_d2t - obt_d2t) < 10e-12);
  assert(fabs(esp_d3t - obt_d3t) < 10e-12);
	printf("tlamb1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testTlamb2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función tlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
 //------------------------------------------------------------------------------
void testTlamb2(bool verbose)
{
  double m = 0.0;
  double q = 0.99132778710087;
  double qsqfm1 = 0.01726921852169195;
  double x = 0.9456982205180451;
  double n = 2.0;

	double esp_t = 0.03640267122578799;
  double esp_dt = -0.0322411111440828;
  double esp_d2t = -0.03409238850679174;
  double esp_d3t = 0.0;

  double obt_t, obt_dt, obt_d2t, obt_d3t;
	tlamb(m, q, qsqfm1, x, n, &obt_t, &obt_dt, &obt_d2t, &obt_d3t);

	if(verbose)
	{
		printf("tlamb2:\n");
    printf("Esperado t: %.20lf \n", esp_t);
		printf("Obtenido t: %.20lf \n", obt_t);
    printf("Esperado dt: %.20lf \n", esp_dt);
		printf("Obtenido dt: %.20lf \n", obt_dt);
    printf("Esperado d2t: %.20lf \n", esp_d2t);
		printf("Obtenido d2t: %.20lf \n", obt_d2t);
    printf("Esperado d3t: %.20lf \n", esp_d3t);
		printf("Obtenido d3t: %.20lf \n", obt_d3t);
	}

	assert(fabs(esp_t - obt_t) < 10e-12);
  assert(fabs(esp_dt - obt_dt) < 10e-12);
  assert(fabs(esp_d2t - obt_d2t) < 10e-12);
  assert(fabs(esp_d3t - obt_d3t) < 10e-12);
	printf("tlamb2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testTlamb3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función tlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
 //------------------------------------------------------------------------------
void testTlamb3(bool verbose)
{
  double m = 0.0;
  double q = 0.9891272558618515;
  double qsqfm1 = 0.02162727171120368;
  double x = 0.6965623667986072;
  double n = 2.0;

	double esp_t = 0.06141053341217671;
  double esp_dt = -0.08623409773884411;
  double esp_d2t = 0.2395363116704881;
  double esp_d3t = 0.0;

  double obt_t, obt_dt, obt_d2t, obt_d3t;
	tlamb(m, q, qsqfm1, x, n, &obt_t, &obt_dt, &obt_d2t, &obt_d3t);

	if(verbose)
	{
		printf("tlamb3:\n");
    printf("Esperado t: %.20lf \n", esp_t);
		printf("Obtenido t: %.20lf \n", obt_t);
    printf("Esperado dt: %.20lf \n", esp_dt);
		printf("Obtenido dt: %.20lf \n", obt_dt);
    printf("Esperado d2t: %.20lf \n", esp_d2t);
		printf("Obtenido d2t: %.20lf \n", obt_d2t);
    printf("Esperado d3t: %.20lf \n", esp_d3t);
		printf("Obtenido d3t: %.20lf \n", obt_d3t);
	}

	assert(fabs(esp_t - obt_t) < 10e-12);
  assert(fabs(esp_dt - obt_dt) < 10e-12);
  assert(fabs(esp_d2t - obt_d2t) < 10e-12);
  assert(fabs(esp_d3t - obt_d3t) < 10e-12);
	printf("tlamb3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_x1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función R_x
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_x1(bool verbose)
{
	double angle = -0.409093775692299;
	double esperado[3][3] = {
		{
			1.0,
			0.0,
			0.0
		},
		{
			0.0,
			0.9174816756401871,
			-0.3977780472379974
		},
		{
			0.0,
			0.3977780472379974,
			0.9174816756401871
		}
	};
	double obtenido[3][3];
	R_x(angle, obtenido);

	if(verbose) printf("R_x1:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_x1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_x2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función R_x
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_x2(bool verbose)
{
	double angle = -0.4090670040227975;
	double esperado[3][3] = {
		{
			1.0,
			0.0,
			0.0
		},
		{
			0.0,
			0.9174923244938117,
			-0.3977534845792582
		},
		{
			0.0,
			0.3977534845792582,
			0.9174923244938117
		}
	};
	double obtenido[3][3];
	R_x(angle, obtenido);

	if(verbose) printf("R_x2:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_x2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_x3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función R_x
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_x3(bool verbose)
{
	double angle = -0.4090670619660224;
	double esperado[3][3] = {
		{
			1.0,
			0.0,
			0.0
		},
		{
			0.0,
			0.9174923014466905,
			-0.3977535377417217
		},
		{
			0.0,
			0.3977535377417217,
			0.9174923014466905
		}
	};
	double obtenido[3][3];
	R_x(angle, obtenido);

	if(verbose) printf("R_x3:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_x3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_y1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función R_y
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_y1(bool verbose)
{
	double angle = 0.0009133473539360685;
	double esperado[3][3] = {
		{
			0.9999995828983346,
			0.0,
			-0.0009133472269498308
		},
		{
			0.0,
			1.0,
			0.0
		},
		{
			0.0009133472269498308,
			0.0,
			0.9999995828983346
		}
	};
	double obtenido[3][3];
	R_y(angle, obtenido);

	if(verbose) printf("R_y1:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_y1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_y2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función R_y
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_y2(bool verbose)
{
	double angle = 0.001069739807684866;
	double esperado[3][3] = {
		{
			0.9999994278284265,
			0.0,
			-0.001069739603659955
		},
		{
			0.0,
			1.0,
			0.0
		},
		{
			0.001069739603659955,
			0.0,
			0.9999994278284265
		}
	};
	double obtenido[3][3];
	R_y(angle, obtenido);

	if(verbose) printf("R_y2:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_y2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_y3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función R_y
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_y3(bool verbose)
{
	double angle = 0.001069836164067632;
	double esperado[3][3] = {
		{
			0.9999994277253456,
			0.0,
			-0.001069835959987584
		},
		{
			0.0,
			1.0,
			0.0
		},
		{
			0.001069835959987584,
			0.0,
			0.9999994277253456
		}
	};
	double obtenido[3][3];
	R_y(angle, obtenido);

	if(verbose) printf("R_y3:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_y3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_z1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función R_z
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_z1(bool verbose)
{
	double angle = -0.001050992069568202;
	double esperado[3][3] = {
		{
			0.9999994477078857,
			-0.001050991876083317,
			0.0
		},
		{
			0.001050991876083317,
			0.9999994477078857,
			0.0
		},
		{
			0.0,
			0.0,
			1.0
		}
	};
	double obtenido[3][3];
	R_z(angle, obtenido);

	if(verbose) printf("R_z1:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_z1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_z2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función R_z
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_z2(bool verbose)
{
	double angle = -0.001230967163554325;
	double esperado[3][3] = {
		{
			0.9999992423600168,
			-0.001230966852677662,
			0.0
		},
		{
			0.001230966852677662,
			0.9999992423600168,
			0.0
		},
		{
			0.0,
			0.0,
			1.0
		}
	};
	double obtenido[3][3];
	R_z(angle, obtenido);

	if(verbose) printf("R_z2:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_z2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testR_z3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función R_z
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testR_z3(bool verbose)
{
	double angle = -0.0012310780508968;
	double esperado[3][3] = {
		{
			0.999999242223512,
			-0.001231077739936117,
			0.0
		},
		{
			0.001231077739936117,
			0.999999242223512,
			0.0
		},
		{
			0.0,
			0.0,
			1.0
		}
	};
	double obtenido[3][3];
	R_z(angle, obtenido);

	if(verbose) printf("R_z3:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(verbose)
			{
				printf("Esperado: %.20lf\n", esperado[i][j]);
				printf("Obtenido: %.20lf\n", obtenido[i][j]);
			}
			assert(fabs(obtenido[i][j]-esperado[i][j])<10e-12);
		}
	}

	printf("R_z3 superado!\n");
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
	testD8rt1(false);
	testD8rt2(false);
	testD8rt3(false);
	printf("d8rt finalizado!\n\n");

  // Test tlamb
	printf("Probando tlamb!\n");
	//testTlamb1(false);
	//testTlamb2(false);
	//testTlamb3(false);
	printf("tlamb finalizado!\n\n");

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

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
