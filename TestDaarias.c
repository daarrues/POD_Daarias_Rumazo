#include "angl.h"
#include "gibbs.h"
#include "lambert_gooding.h"
#include "Mjday.h"
#include "newtonnu.h"
#include "Position.h"
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


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
//  void testXlamb1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función xlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testXlamb1(bool verbose)
{
  double m = 0.0;
  double q = 0.8046115644662321;
  double qsqfm1 = 0.3526002303272028;
  double tin = 0.9134381609833159;

	double esp_n = 1.0;
  double esp_x = 0.6256383152410641;
  double esp_xpl = 0.0;

  double obt_n, obt_x, obt_xpl;
	xlamb(m, q, qsqfm1, tin, &obt_n, &obt_x, &obt_xpl);

	if(verbose)
	{
		printf("xlamb1:\n");
    printf("Esperado n: %.20lf \n", esp_n);
		printf("Obtenido n: %.20lf \n", obt_n);
    printf("Esperado x: %.20lf \n", esp_x);
		printf("Obtenido x: %.20lf \n", obt_x);
    printf("Esperado xpl: %.20lf \n", esp_xpl);
		printf("Obtenido xpl: %.20lf \n", obt_xpl);
	}

	assert(fabs(esp_n - obt_n) < 10e-12);
  assert(fabs(esp_x - obt_x) < 10e-12);
  assert(fabs(esp_xpl - obt_xpl) < 10e-12);
	printf("xlamb1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testXlamb2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función xlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testXlamb2(bool verbose)
{
  double m = 1.0;
  double q = 0.99132778710087;
  double qsqfm1 = 0.01726921852169195;
  double tin = 0.0637845039365542;

	double esp_n = 0.0;
  double esp_x = 0.0;
  double esp_xpl = 0.0;

  double obt_n, obt_x, obt_xpl;
	xlamb(m, q, qsqfm1, tin, &obt_n, &obt_x, &obt_xpl);

	if(verbose)
	{
		printf("xlamb2:\n");
    printf("Esperado n: %.20lf \n", esp_n);
		printf("Obtenido n: %.20lf \n", obt_n);
    printf("Esperado x: %.20lf \n", esp_x);
		printf("Obtenido x: %.20lf \n", obt_x);
    printf("Esperado xpl: %.20lf \n", esp_xpl);
		printf("Obtenido xpl: %.20lf \n", obt_xpl);
	}

	assert(fabs(esp_n - obt_n) < 10e-12);
  assert(fabs(esp_x - obt_x) < 10e-12);
  assert(fabs(esp_xpl - obt_xpl) < 10e-12);
	printf("xlamb2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testXlamb3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función xlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testXlamb3(bool verbose)
{
  double m = 0.0;
  double q = 0.9669569195016101;
  double qsqfm1 = 0.06499431582795664;
  double tin = 0.1741579527629812;

	double esp_n = 1.0;
  double esp_x = 0.7236378712038933;
  double esp_xpl = 0.0;

  double obt_n, obt_x, obt_xpl;
	xlamb(m, q, qsqfm1, tin, &obt_n, &obt_x, &obt_xpl);

	if(verbose)
	{
		printf("xlamb3:\n");
    printf("Esperado n: %.20lf \n", esp_n);
		printf("Obtenido n: %.20lf \n", obt_n);
    printf("Esperado x: %.20lf \n", esp_x);
		printf("Obtenido x: %.20lf \n", obt_x);
    printf("Esperado xpl: %.20lf \n", esp_xpl);
		printf("Obtenido xpl: %.20lf \n", obt_xpl);
	}

	assert(fabs(esp_n - obt_n) < 10e-12);
  assert(fabs(esp_x - obt_x) < 10e-12);
  assert(fabs(esp_xpl - obt_xpl) < 10e-12);
	printf("xlamb3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testVlamb1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función vlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testVlamb1(bool verbose)
{
  double gm = 398600441800000;
  double r1 = 9163781.535461562;
  double r2 = 9158926.303943245;
  double th = 0.431406120912815;
	double tdelt = 600.0000044703484;

	double esp_n = 1.0;
  double esp_vri[2] = {
		-4.246341652783283,
		0.0
		};
	double esp_vti[2] = {
		6585.897930820493,
		0.0
		};
	double esp_vrf[2] = {
		-11.68727565166062,
		0.0
		};
	double esp_vtf[2] = {
		6589.38917620767,
		0.0
		};

  double obt_n, obt_vri[2], obt_vti[2], obt_vrf[2], obt_vtf[2];
	vlamb(gm, r1, r2, th, tdelt, &obt_n, obt_vri, obt_vti, obt_vrf, obt_vtf);

	if(verbose){
		printf("vlamb1:\n");

		// N
		printf("Esperado n: %.20lf \n", esp_n);
		printf("Obtenido n: %.20lf \n", obt_n);
	}
	assert(fabs(esp_n - obt_n) < 10e-12);

	// VRI
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vri[i]);
			printf("Obtenido: %.20lf\n", obt_vri[i]);
		}
		assert(fabs(obt_vri[i]-esp_vri[i])<10e-12);
	}

	// VTI
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vti[i]);
			printf("Obtenido: %.20lf\n", obt_vti[i]);
		}
		assert(fabs(obt_vti[i]-esp_vti[i])<10e-12);
	}

	// VRF
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vrf[i]);
			printf("Obtenido: %.20lf\n", obt_vrf[i]);
		}
		assert(fabs(obt_vrf[i]-esp_vrf[i])<10e-12);
	}

	// VTF
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vtf[i]);
			printf("Obtenido: %.20lf\n", obt_vtf[i]);
		}
		assert(fabs(obt_vtf[i]-esp_vtf[i])<10e-12);
	}

	printf("vlamb1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testVlamb2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función vlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testVlamb2(bool verbose)
{
  double gm = 398600441800000;
  double r1 = 9163781.535461562;
  double r2 = 9158926.303943245;
  double th = 6.714591428092401;
	double tdelt = 600.0000044703484;

	double esp_n = 0.0;
  double esp_vri[2] = {
		0.0,
		0.0
		};
	double esp_vti[2] = {
		0.0,
		0.0
		};
	double esp_vrf[2] = {
		0.0,
		0.0
		};
	double esp_vtf[2] = {
		0.0,
		0.0
		};

  double obt_n, obt_vri[2], obt_vti[2], obt_vrf[2], obt_vtf[2];
	vlamb(gm, r1, r2, th, tdelt, &obt_n, obt_vri, obt_vti, obt_vrf, obt_vtf);

	if(verbose)
	{
		printf("vlamb2:\n");

		// N
		printf("Esperado n: %.20lf \n", esp_n);
		printf("Obtenido n: %.20lf \n", obt_n);
	}
	assert(fabs(esp_n - obt_n) < 10e-12);

	// VRI
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vri[i]);
			printf("Obtenido: %.20lf\n", obt_vri[i]);
		}
		assert(fabs(obt_vri[i]-esp_vri[i])<10e-12);
	}

	// VTI
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vti[i]);
			printf("Obtenido: %.20lf\n", obt_vti[i]);
		}
		assert(fabs(obt_vti[i]-esp_vti[i])<10e-12);
	}

	// VRF
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vrf[i]);
			printf("Obtenido: %.20lf\n", obt_vrf[i]);
		}
		assert(fabs(obt_vrf[i]-esp_vrf[i])<10e-12);
	}

	// VTF
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vtf[i]);
			printf("Obtenido: %.20lf\n", obt_vtf[i]);
		}
		assert(fabs(obt_vtf[i]-esp_vtf[i])<10e-12);
	}

	printf("vlamb2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testVlamb3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función vlamb
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testVlamb3(bool verbose)
{
  double gm = 398600441800000;
  double r1 = 20473061.83367915451526641846;
  double r2 = 20488505.59583704546093940735;
  double th = 0.06720882293144518627;
	double tdelt = 300.00002235174179077148;

	double esp_n = 1.0;
  double esp_vri[2] = {
		39.82208832992758829050,
		0.0
		};
	double esp_vti[2] = {
		4589.76670044189813779667,
		0.0
		};
	double esp_vrf[2] = {
		63.09171280899065692438,
		0.0
		};
	double esp_vtf[2] = {
		4586.30704034369045984931,
		0.0
		};

  double obt_n, obt_vri[2], obt_vti[2], obt_vrf[2], obt_vtf[2];
	vlamb(gm, r1, r2, th, tdelt, &obt_n, obt_vri, obt_vti, obt_vrf, obt_vtf);

	if(verbose)
	{
		printf("vlamb3:\n");

		// N
		printf("Esperado n: %.20lf \n", esp_n);
		printf("Obtenido n: %.20lf \n", obt_n);
	}
	assert(fabs(esp_n - obt_n) < 10e-12);

	// VRI
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vri[i]);
			printf("Obtenido: %.20lf\n", obt_vri[i]);
		}
		assert(fabs(obt_vri[i]-esp_vri[i])<10e-12);
	}

	// VTI
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vti[i]);
			printf("Obtenido: %.20lf\n", obt_vti[i]);
		}
		assert(fabs(obt_vti[i]-esp_vti[i])<10e-12);
	}

	// VRF
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vrf[i]);
			printf("Obtenido: %.20lf\n", obt_vrf[i]);
		}
		assert(fabs(obt_vrf[i]-esp_vrf[i])<10e-12);
	}

	// VTF
	for(int i = 0; i < 2; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_vtf[i]);
			printf("Obtenido: %.20lf\n", obt_vtf[i]);
		}
		assert(fabs(obt_vtf[i]-esp_vtf[i])<10e-12);
	}

	printf("vlamb3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testLambert_gooding1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función lambert_gooding
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testLambert_gooding1(bool verbose)
{
	double r1[3] = {
		8794276.580984011,
		404708.1949434897,
		2543973.805637163
		};
	double r2[3] = {
		8330586.996205064,
		3762923.082128749,
		572416.9963455135
		};
  double tof = 600.0000044703484;
  double mu = 398600441800000;
  double long_way = 0.0;
  double multi_revs = 1.0;

	double esp_v1[3] = {
		591.4156797181868,
		5838.863650451791,
		-2988.63988326117
		};
	double esp_v2[3] = {
		-2113.653770404029,
		5180.392996064316,
		-3480.83071307979
		};


  double obt_v1[3], obt_v2[3];
	lambert_gooding(r1, r2, tof, mu, long_way, multi_revs, obt_v1, obt_v2);

	if(verbose) printf("lambert_gooding1:\n");

	// V1
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_v1[i]);
			printf("Obtenido: %.20lf\n", obt_v1[i]);
		}
		assert(fabs(obt_v1[i]-esp_v1[i])<10e-12);
	}

	// V2
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_v2[i]);
			printf("Obtenido: %.20lf\n", obt_v2[i]);
		}
		assert(fabs(obt_v2[i]-esp_v2[i])<10e-12);
	}

	printf("lambert_gooding1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testLambert_gooding2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función lambert_gooding
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testLambert_gooding2(bool verbose)
{
	double r1[3] = {
		1976694.85762782581150531769,
		-37017163.92619637399911880493,
		-17995197.21750879287719726562
		};
	double r2[3] = {
		2464157.77993164770305156708,
		-36524015.15769488364458084106,
		-17826894.45484596118330955505
		};
  double tof = 299.99998211860656738281;
  double mu = 398600441800000;
  double long_way = 0.0;
  double multi_revs = 1.0;

	double esp_v1[3] = {
		1493.24026997438159014564,
		1473.89412066475870233262,
		497.96118943721626237675
		};
	double esp_v2[3] = {
		1489.02795936263919429621,
		1543.60206023751629800245,
		531.91650817821096097759
		};


  double obt_v1[3], obt_v2[3];
	lambert_gooding(r1, r2, tof, mu, long_way, multi_revs, obt_v1, obt_v2);

	if(verbose) printf("lambert_gooding2:\n");

	// V1
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_v1[i]);
			printf("Obtenido: %.20lf\n", obt_v1[i]);
		}
		assert(fabs(obt_v1[i]-esp_v1[i])<10e-9); // 12 significativas (al menos)
	}

	// V2
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_v2[i]);
			printf("Obtenido: %.20lf\n", obt_v2[i]);
		}
		assert(fabs(obt_v2[i]-esp_v2[i])<10e-9); // 12 significativas (al menos)
	}

	printf("lambert_gooding2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testLambert_gooding3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función lambert_gooding
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testLambert_gooding3(bool verbose)
{
	double r1[3] = {
		20387627.07175297,
		1865163.696333993,
		-109943.6885558921
		};
	double r2[3] = {
		20435422.35215281,
		1070699.446717979,
		1012905.491433884
		};
  double tof = 300.0000223517418;
  double mu = 398600441800000;
  double long_way = 0.0;
  double multi_revs = 1.0;

	double esp_v1[3] = {
		301.4349823754282,
		-2637.066656808077,
		3744.67095512485
		};
	double esp_v2[3] = {
		17.19646978067542,
		-2657.510276115648,
		3738.386850808526
		};


  double obt_v1[3], obt_v2[3];
	lambert_gooding(r1, r2, tof, mu, long_way, multi_revs, obt_v1, obt_v2);

	if(verbose) printf("lambert_gooding3:\n");

	// V1
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_v1[i]);
			printf("Obtenido: %.20lf\n", obt_v1[i]);
		}
		assert(fabs(obt_v1[i]-esp_v1[i])<10e-9); // 12 significativas (al menos)
	}

	// V2
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf\n", esp_v2[i]);
			printf("Obtenido: %.20lf\n", obt_v2[i]);
		}
		assert(fabs(obt_v2[i]-esp_v2[i])<10e-10); // 12 significativas (al menos)
	}

	printf("lambert_gooding3 superado!\n");
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
//  void testAngl1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función angl
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testAngl1(bool verbose)
{
	double v1[3] = {
		20387627.07175297,
		1865163.696333993,
		-109943.6885558921
		};
	double v2[3] = {
		20435422.35215281,
		1070699.446717979,
		1012905.491433884
		};
	double esperado = 0.06720882293144352;
	double obtenido = angl(v1, v2);

	if(verbose)
	{
		printf("angl1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("angl1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testAngl2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función angl
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testAngl2(bool verbose)
{
	double v1[3] = {
		20435422.35215281,
		1070699.446717979,
		1012905.491433884
		};
	double v2[3] = {
		20398157.06662507,
		271778.6158697309,
		2131538.395420803
		};
	double esperado = 0.06708423348977346;
	double obtenido = angl(v1, v2);

	if(verbose)
	{
		printf("angl2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("angl2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testAngl3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función angl
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testAngl3(bool verbose)
{
	double v1[3] = {
		76378095813.3326,
		6694495484.997974,
		0.0
		};
	double v2[3] = {
		0.08111985681453468,
		0.01287140290519155,
		-0.008099964956862108
		};
	double esperado = 0.1205731168766288;
	double obtenido = angl(v1, v2);

	if(verbose)
	{
		printf("angl3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("angl3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testMjday1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función Mjday
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testMjday1(bool verbose)
{
	double year = 2009.0;
	double month = 5.0;
	double day = 26.0;
	double hour = 16.0;
	double min = 0.0;
	double sec = 20.475;
	double esperado = 54977.66690364573;
	double obtenido = Mjday(year, month, day, hour, min, sec);

	if(verbose)
	{
		printf("Mjday1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("Mjday1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testMjday2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función Mjday
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testMjday2(bool verbose)
{
	double year = 2011.0;
	double month = 1.0;
	double day = 4.0;
	double hour = 13.0;
	double min = 0.0;
	double sec = 46.5;
	double esperado = 55565.54220486106;
	double obtenido = Mjday(year, month, day, hour, min, sec);

	if(verbose)
	{
		printf("Mjday2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("Mjday2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testMjday3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función Mjday
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testMjday3(bool verbose)
{
	double year = 2011.0;
	double month = 1.0;
	double day = 4.0;
	double hour = 21.0;
	double min = 42.0;
	double sec = 20.796;
	double esperado = 55565.90440736106;
	double obtenido = Mjday(year, month, day, hour, min, sec);

	if(verbose)
	{
		printf("Mjday3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-12);
	printf("Mjday3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNewtonnu1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función newtonnu
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNewtonnu1(bool verbose)
{
	double ecc = 0.08253310617358817;
	double nu = 0.1812003110535373;
	double esp_e0 = 0.1668839377415842;
	double esp_m = 0.153174331355188;

	double obt_e0, obt_m;
	newtonnu(ecc, nu, &obt_e0, &obt_m);

	if(verbose)
	{
		printf("newtonnu1:\n");
		printf("Esperado e0: %.20lf \n", esp_e0);
		printf("Obtenido e0: %.20lf \n", obt_e0);
		printf("Esperado m: %.20lf \n", esp_m);
		printf("Obtenido m: %.20lf \n", obt_m);
	}

	assert(fabs(esp_e0 - obt_e0) < 10e-12);
	assert(fabs(esp_m - obt_m) < 10e-12);
	printf("newtonnu1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNewtonnu2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función newtonnu
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNewtonnu2(bool verbose)
{
	double ecc = 0.08669431216172874;
	double nu = 0.1713058522671266;
	double esp_e0 = 0.157107124092798;
	double esp_m = 0.1435427917466885;

	double obt_e0, obt_m;
	newtonnu(ecc, nu, &obt_e0, &obt_m);

	if(verbose)
	{
		printf("newtonnu2:\n");
		printf("Esperado e0: %.20lf \n", esp_e0);
		printf("Obtenido e0: %.20lf \n", obt_e0);
		printf("Esperado m: %.20lf \n", esp_m);
		printf("Obtenido m: %.20lf \n", obt_m);
	}

	assert(fabs(esp_e0 - obt_e0) < 10e-12);
	assert(fabs(esp_m - obt_m) < 10e-12);
	printf("newtonnu2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testNewtonnu3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función newtonnu
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testNewtonnu3(bool verbose)
{
	double ecc = 0.07912910777801811;
	double nu = 0.1902672374573671;
	double esp_e0 = 0.1758403692885834;
	double esp_m = 0.1619978705540214;

	double obt_e0, obt_m;
	newtonnu(ecc, nu, &obt_e0, &obt_m);

	if(verbose)
	{
		printf("newtonnu3:\n");
		printf("Esperado e0: %.20lf \n", esp_e0);
		printf("Obtenido e0: %.20lf \n", obt_e0);
		printf("Esperado m: %.20lf \n", esp_m);
		printf("Obtenido m: %.20lf \n", obt_m);
	}

	assert(fabs(esp_e0 - obt_e0) < 10e-12);
	assert(fabs(esp_m - obt_m) < 10e-12);
	printf("newtonnu1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPosition1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función Position
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPosition1(bool verbose)
{
	double lon = -2.117969613665733;
	double lat = 0.6830532777909771;
	double h = 99.81638;

	double esperado[3] = {
		-2577383.639573138,
		-4230610.246798722,
		4004108.332058704
		};

	double obtenido[3];
	Position(lon, lat, h, obtenido);

	if(verbose) printf("Position1:\n");
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf \n", esperado[i]);
			printf("Obtenido: %.20lf \n", obtenido[i]);
		}
		assert(fabs(esperado[i] - obtenido[i]) < 10e-9);
		// 12 significativas (al menos)
	}

	printf("Position1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPosition2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función Position
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPosition2(bool verbose)
{
	double lon = -1.504723397302147;
	double lat = 0.5335890402367144;
	double h = 0.0;

	double esperado[3] = {
		362889.5147507534,
		-5484262.361013475,
		3225167.728477614
		};

	double obtenido[3];
	Position(lon, lat, h, obtenido);

	if(verbose) printf("Position2:\n");
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf \n", esperado[i]);
			printf("Obtenido: %.20lf \n", obtenido[i]);
		}
		assert(fabs(esperado[i] - obtenido[i]) < 10e-9);
		// 12 significativas (al menos)
	}

	printf("Position2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testPosition3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función Position
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testPosition3(bool verbose)
{
	double lon = -1.504723397302147;
	double lat = 0.5335890402367144;
	double h = 0.0;

	double esperado[3] = {
		362889.5147507534,
		-5484262.361013475,
		3225167.728477614
		};

	double obtenido[3];
	Position(lon, lat, h, obtenido);

	if(verbose) printf("Position3:\n");
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf \n", esperado[i]);
			printf("Obtenido: %.20lf \n", obtenido[i]);
		}
		assert(fabs(esperado[i] - obtenido[i]) < 10e-9);
		// 12 significativas (al menos)
	}

	printf("Position3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGibbs1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función gibbs
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGibbs1(bool verbose)
{
	double r1[3] = {
		20387627.07175297,
		1865163.696333993,
		-109943.6885558921
		};
	double r2[3] = {
		20435422.35215281,
		1070699.446717979,
		1012905.491433884
		};
	double r3[3] = {
		20398157.06662507,
		271778.6158697309,
		2131538.395420803
		};

	double esp_v2[3] = {
		17.44484600828158,
		-2659.686950246493,
		3741.477704717453
		};
	double esp_theta = 0.06720882293144352;
	double esp_theta1 = 0.06708423348977346;
	double esp_copa = 1.737759242059767e-15err;
	char esp_error[12] = "          ok";

	double obt_v2[3], obt_theta, obt_theta1, obt_copa;
	char obt_error[12];
	gibbs(r1, r2, r3, obt_v2, &obt_theta, &obt_theta1, &obt_copa, obt_error);

	if(verbose) printf("gibbs1:\n");
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf \n", esp_v2[i]);
			printf("Obtenido: %.20lf \n", obt_v2[i]);
		}
		assert(fabs(esp_v2[i] - obt_v2[i]) < 10e-12);
	}
	if(verbose){
		printf("Esperado: %.20lf \n", esp_theta[i]);
		printf("Obtenido: %.20lf \n", obt_theta[i]);
		printf("Esperado: %.20lf \n", esp_theta1[i]);
		printf("Obtenido: %.20lf \n", obt_theta1[i]);
		printf("Esperado: %.20lf \n", esp_copa[i]);
		printf("Obtenido: %.20lf \n", obt_copa[i]);
		printf("Esperado: %s \n", esp_error);
		printf("Obtenido: %s \n", obt_error);
	}
	assert(fabs(esp_theta[i] - obt_theta[i]) < 10e-12);
	assert(fabs(esp_theta1[i] - obt_theta1[i]) < 10e-12);
	assert(fabs(esp_copa[i] - obt_copa[i]) < 10e-12);
	if(strcmp(esp_error, obt_error) != 0) assert(2 < 1);

	printf("gibbs1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGibbs2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función gibbs
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGibbs2(bool verbose)
{
	double r1[3] = {
		20408482.65293036,
		1869905.834429736,
		-114536.6708983118
		};
	double r2[3] = {
		20456329.59101998,
		1074191.366837643,
		1009857.021678711
		};
	double r3[3] = {
		20418881.02782273,
		273996.4780377855,
		2130041.987544681
		};

	double esp_v2[3] = {
		17.1865619482012,
		-2657.52837347195,
		3737.684787287446
		};
	double esp_theta = 0.06723676688978014;
	double esp_theta1 = 0.06711481197035603;
	double esp_copa = 1.992329912159363e-15;
	char esp_error[12] = "          ok";

	double obt_v2[3], obt_theta, obt_theta1, obt_copa;
	char obt_error[12];
	gibbs(r1, r2, r3, obt_v2, &obt_theta, &obt_theta1, &obt_copa, obt_error);

	if(verbose) printf("gibbs2:\n");
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf \n", esp_v2[i]);
			printf("Obtenido: %.20lf \n", obt_v2[i]);
		}
		assert(fabs(esp_v2[i] - obt_v2[i]) < 10e-12);
	}
	if(verbose){
		printf("Esperado: %.20lf \n", esp_theta[i]);
		printf("Obtenido: %.20lf \n", obt_theta[i]);
		printf("Esperado: %.20lf \n", esp_theta1[i]);
		printf("Obtenido: %.20lf \n", obt_theta1[i]);
		printf("Esperado: %.20lf \n", esp_copa[i]);
		printf("Obtenido: %.20lf \n", obt_copa[i]);
		printf("Esperado: %s \n", esp_error);
		printf("Obtenido: %s \n", obt_error);
	}
	assert(fabs(esp_theta[i] - obt_theta[i]) < 10e-12);
	assert(fabs(esp_theta1[i] - obt_theta1[i]) < 10e-12);
	assert(fabs(esp_copa[i] - obt_copa[i]) < 10e-12);
	if(strcmp(esp_error, obt_error) != 0) assert(2 < 1);

	printf("gibbs2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testGibbs3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función gibbs
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testGibbs3(bool verbose)
{
	double r1[3] = {
		20370508.34273631,
		1861271.242974382,
		-106173.665585571
		};
	double r2[3] = {
		20418280.37423369,
		1067836.399236342,
		1015404.951145533
		};
	double r3[3] = {
		20381175.29424456,
		269961.2398309496,
		2132764.59236734
		};

	double esp_v2[3] = {
		17.70552956104248,
		-2661.30974048117,
		3744.38790480311
		};
	double esp_theta = 0.06718562571898608;
	double esp_theta1 = 0.06705901381744461;
	double esp_copa = 1.634109514370152e-15;
	char esp_error[12] = "          ok";

	double obt_v2[3], obt_theta, obt_theta1, obt_copa;
	char obt_error[12];
	gibbs(r1, r2, r3, obt_v2, &obt_theta, &obt_theta1, &obt_copa, obt_error);

	if(verbose) printf("gibbs3:\n");
	for(int i = 0; i < 3; i++)
	{
		if(verbose)
		{
			printf("Esperado: %.20lf \n", esp_v2[i]);
			printf("Obtenido: %.20lf \n", obt_v2[i]);
		}
		assert(fabs(esp_v2[i] - obt_v2[i]) < 10e-12);
	}
	if(verbose){
		printf("Esperado: %.20lf \n", esp_theta[i]);
		printf("Obtenido: %.20lf \n", obt_theta[i]);
		printf("Esperado: %.20lf \n", esp_theta1[i]);
		printf("Obtenido: %.20lf \n", obt_theta1[i]);
		printf("Esperado: %.20lf \n", esp_copa[i]);
		printf("Obtenido: %.20lf \n", obt_copa[i]);
		printf("Esperado: %s \n", esp_error);
		printf("Obtenido: %s \n", obt_error);
	}
	assert(fabs(esp_theta[i] - obt_theta[i]) < 10e-12);
	assert(fabs(esp_theta1[i] - obt_theta1[i]) < 10e-12);
	assert(fabs(esp_copa[i] - obt_copa[i]) < 10e-12);
	if(strcmp(esp_error, obt_error) != 0) assert(2 < 1);

	printf("gibbs3 superado!\n");
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

	// Test angl
	printf("Probando angl!\n");
	testAngl1(false);
	testAngl2(false);
	testAngl3(false);
	printf("angl finalizado!\n\n");

	// Test Mjday
	printf("Probando Mjday!\n");
	testMjday1(false);
	testMjday2(false);
	testMjday3(false);
	printf("Mjday finalizado!\n\n");

	// Test newtonnu
	printf("Probando newtonnu!\n");
	testNewtonnu1(false);
	testNewtonnu2(false);
	testNewtonnu3(false);
	printf("newtonnu finalizado!\n\n");

	// Test Position
	printf("Probando Position!\n");
	testPosition1(false);
	testPosition2(false);
	testPosition3(false);
	printf("Position finalizado!\n\n");

	// Test gibbs
	printf("Probando Position!\n");
	testGibbs1(false);
	testGibbs2(false);
	testGibbs3(false);
	printf("Position finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
