#include "MatLabUtils.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

typedef int bool;
#define true 1
#define false 0
#define POL_DEG 15

void testNorm1(bool verbose)
{
	double v[] = {
		4950990.33826460,
		256563.116260381,
		3999465.34658133
		};
	double esperado = 6.36976082916144930000e+06;
	double obtenido = norm(v);

	if(verbose)
	{
		printf("Norm1:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-6);
	printf("Norm1 superado!\n");
}

void testNorm2(bool verbose)
{
	double v[] = {
		4935037.85913036,
		472703.320202615,
		3999475.70573182,
		};
	double esperado = 6.36976082916145030000e+06;
	double obtenido = norm(v);

	if(verbose)
	{
		printf("Norm2:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-6);
	printf("Norm2 superado!\n");
}

void testNorm3(bool verbose)
{
	double v[] = {
		4909646.95198536,
		687938.936915757,
		3999494.94894739,
		};
	double esperado = 6.36976082916144840000e+06;
	double obtenido = norm(v);

	if(verbose)
	{
		printf("Norm3:\n");
		printf("Esperado: %.20lf \n", esperado);
		printf("Obtenido: %.20lf \n", obtenido);
	}

	assert(fabs(esperado - obtenido) < 10e-6);
	printf("Norm3 superado!\n");
}

void testDot1(bool verbose)
{
	double v1[] = {
		6.33886095165153710000e-01,
		-7.73340379778256310000e-01,
		1.15358294325144310000e-02
		};

	double v2[] = {
		4.95099033826459660000e+06,
		2.56563116260380720000e+05,
		3.99946534658133240000e+06
		};

	double expected = 2.98609046501646750000e+06;

	double res = dot(v1,v2);

	if(verbose)
	{
		printf("Dot1:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}

	assert(fabs(res-expected)<10e-12);
	printf("Dot1 superado!\n");
}

void testDot2(bool verbose)
{
	double v1[] = {
		9.35539825569648990000e-01,
		-1.64830818224054180000e-02,
		-3.52836424971610050000e-01
		};

	double v2[] = {
		4.93503785913035740000e+06,
		4.72703320202615400000e+05,
		3.99947570573182080000e+06
		};

	double expected = 3.19797214063458420000e+06;

	double res = dot(v1,v2);

	if(verbose)
	{
		printf("Dot2:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}

	assert(fabs(res-expected)<10e-12);
	printf("Dot2 superado!\n");
}

void testDot3(bool verbose)
{
	double v1[] = {
		2.04354223521544110000e+07,
		1.07069944671824550000e+06,
		1.01290549143365120000e+06
		};

	double v2[] = {
		1.71964697862374220000e+01,
		-2.65751027611478370000e+03,
		3.73838685080780670000e+03
		};

	double expected = 1.29265491105025480000e+09;

	double res = dot(v1,v2);

	if(verbose)
	{
		printf("Dot3:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}

	assert(fabs(res-expected)<10e-12);
	printf("Dot3 superado!\n");
}

void testCross1(bool verbose)
{
	double v1[] = {
		20435422.352153745,
		1070699.4467181342,
		1012905.4914337485
	};

	double v2[] = {
		17.196469780754562,
		-2657.5102761159114,
		3738.386850808869
	};

	double vResult[3];

	double esperado[] = {
		6694495484.9988279,
		-76378095813.343018,
		-54325757148.297119
	};

	cross(v1, v2, vResult);

	if(verbose)
	{
		printf("Cross1:\n");
		printf("Esperado_x = %.20lf\n", esperado[0]);
		printf("Obtenido_x = %.20lf\n", vResult[0]);
		printf("Esperado_y = %.20lf\n", esperado[1]);
		printf("Obtenido_y = %.20lf\n", vResult[1]);
		printf("Esperado_z = %.20lf\n", esperado[2]);
		printf("Obtenido_z = %.20lf\n", vResult[2]);
	}

	assert(fabs(vResult[0]-esperado[0])<10e-12);
	assert(fabs(vResult[1]-esperado[1])<10e-12);
	assert(fabs(vResult[2]-esperado[2])<10e-12);
	printf("Cross1 superado!\n");
}

void testCross2(bool verbose)
{
	double v1[] = {
		20456329.59102045000000000000,
		1074191.36683772060000000000,
		1009857.02167864280000000000
	};

	double v2[] = {
		17.65947568968910300000,
		-2661.67627260588730000000,
		3743.57355524540890000000
	};

	double vResult[3];

	double esperado[] = {
		6709226867.49310400000000000000,
		-76561940948.80389400000000000000,
		-54467096753.35356900000000000000
	};

	cross(v1, v2, vResult);

	if(verbose)
	{
		printf("Cross2:\n");
		printf("Esperado_x = %.20lf\n", esperado[0]);
		printf("Obtenido_x = %.20lf\n", vResult[0]);
		printf("Esperado_y = %.20lf\n", esperado[1]);
		printf("Obtenido_y = %.20lf\n", vResult[1]);
		printf("Esperado_z = %.20lf\n", esperado[2]);
		printf("Obtenido_z = %.20lf\n", vResult[2]);
	}

	assert(fabs(vResult[0]-esperado[0])<10e-12);
	assert(fabs(vResult[1]-esperado[1])<10e-12);
	assert(fabs(vResult[2]-esperado[2])<10e-12);
	printf("Cross2 superado!\n");
}

void testCross3(bool verbose)
{
	double v1[] = {
		20418280.37423648300000000000,
		1067836.39923680880000000000,
		1015404.95114512560000000000
	};

	double v2[] = {
		16.87979502296278200000,
		-2654.08002932726280000000,
		3734.12004615458320000000
	};

	double vResult[3];

	double esperado[] = {
		6682395306.91799930000000000000,
		-76227170226.00053400000000000000,
		-54209775034.00301400000000000000
	};

	cross(v1, v2, vResult);

	if(verbose)
	{
		printf("Cross3:\n");
		printf("Esperado_x = %.20lf\n", esperado[0]);
		printf("Obtenido_x = %.20lf\n", vResult[0]);
		printf("Esperado_y = %.20lf\n", esperado[1]);
		printf("Obtenido_y = %.20lf\n", vResult[1]);
		printf("Esperado_z = %.20lf\n", esperado[2]);
		printf("Obtenido_z = %.20lf\n", vResult[2]);
	}

	assert(fabs(vResult[0]-esperado[0])<10e-12);
	assert(fabs(vResult[1]-esperado[1])<10e-12);
	assert(fabs(vResult[2]-esperado[2])<10e-12);
	printf("Cross3 superado!\n");
}

// Aquí va TestLength

void testZeros(bool verbose)
{
	double v[3];
	zeros(v);

	if(verbose)
	{
		printf("Zeros:\n");
		printf("v[0] = %.20lf\n", v[0]);
		printf("v[1] = %.20lf\n", v[1]);
		printf("v[2] = %.20lf\n", v[2]);
	}

	assert(fabs(v[0]) < 10e-12);
	assert(fabs(v[1]) < 10e-12);
	assert(fabs(v[2]) < 10e-12);
	printf("Zeros superado!\n");
}

void testProdMatr1(bool verbose)
{
// Matriz 1
	double m1[3][3] = {
		{
			0.9999999978959845,
			-5.951700510045733e-05,
			-2.580227134293989e-05
		},
		{
			5.951642956254065e-05,
			0.9999999979801206,
			-2.230590149845583e-05
		},
		{
			2.580359887127566e-05,
			2.230436579244361e-05,
			0.9999999994183447
		}
	};

	// Matriz 2
	double m2[3][3] = {
		{
			0.9999973738023294,
			-0.002101948193683347,
			-0.0009133467225153597
		},
		{
			0.002101948193669182,
			0.9999977909039949,
			-9.599205155674982e-07
		},
		{
			0.0009133467225479574,
			-9.598894989583521e-07,
			0.9999995828983346
		}
	};
// Matriz esperada
	double mE[3][3] = {
		{
			0.9999972230302382,
			-0.002161465038115165,
			-0.0009391489240428393
		},
		{
			0.002161444089662105,
			0.9999976638050796,
			-2.33201718441325e-05
		},
		{
			0.0009391971357440173,
			2.12901891935202e-05,
			0.9999995587276367
		}
	};
// Matriz resultado
	double mR[3][3];
prodMatr(m1, m2, mR);

	if(verbose) printf("ProdMatr1:\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
if(verbose)
			{
				printf("Exp.:%.20lf\n", mE[i][j]);
				printf("Obt.:%.20lf\n", mR[i][j]);
			}
			assert(fabs(mR[i][j]-mE[i][j])<10e-12);
		}
	}
	printf("ProdMatr1 superado!\n");
}

void testTrans(bool verbose)
{
	double m[3][3];
	double t[3][3];

	m[0][0] = 0.307371022581072;
	m[0][1] = 0.951589691503099;
	m[0][2] = -0.000336901326720985;
	m[1][0] = -0.951589111300356;
	m[1][1] = 0.307371207205992;
	m[1][2] = 0.00105082602288526;
	m[2][0] = 0.00110350897844434;
	m[2][1] = -2.40183511869764e-06;
	m[2][2] = 0.999999391130897;

	trans(m, t);

	if(verbose)
	{
		printf("Trans:\n");
		printf("m[0][0] = %.20lf ### t[0][0] = %.20lf\n", m[0][0], t[0][0]);
		printf("m[0][1] = %.20lf ### t[1][0] = %.20lf\n", m[0][1], t[1][0]);
		printf("m[0][2] = %.20lf ### t[2][0] = %.20lf\n", m[0][2], t[2][0]);
		printf("m[1][0] = %.20lf ### t[0][1] = %.20lf\n", m[1][0], t[0][1]);
		printf("m[1][1] = %.20lf ### t[1][1] = %.20lf\n", m[1][1], t[1][1]);
		printf("m[1][2] = %.20lf ### t[2][1] = %.20lf\n", m[1][2], t[2][1]);
		printf("m[2][0] = %.20lf ### t[0][2] = %.20lf\n", m[2][0], t[0][2]);
		printf("m[2][1] = %.20lf ### t[1][2] = %.20lf\n", m[2][1], t[1][2]);
		printf("m[2][2] = %.20lf ### t[2][2] = %.20lf\n", m[2][2], t[2][2]);
	}

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			assert(fabs(m[i][j] - t[j][i]) < 10e-12);
		}
	}
	printf("Trans superado!\n");
}

// Aquí va TestDet

void testRoots(bool verbose)
{
	double p[POL_DEG + 1];
	p[0] = 1.0;
	p[1] = 0.0;
	p[2] = -73120740632113.531;
	p[3] = 0.0;
	p[4] = 0.0;
	p[5] = -1587936795677001700000000000000000000.0;
	p[6] = 0.0;
	p[7] = 0.0;
	p[8] = -11985384853690914000000000000000000000000000000000000000000.0;
	p[9] = 0.0;
	p[10] = 0.0;
	p[11] = 0.0;
	p[12] = 0.0;
	p[13] = 0.0;
	p[14] = 0.0;
	p[15] = 0.0;

	int nRootsEsperado = 9;
	double esperado[9];
	esperado[0] = 0.0;
	esperado[1] = 0.0;
	esperado[2] = 0.0;
	esperado[3] = 0.0;
	esperado[4] = 0.0;
	esperado[5] = 0.0;
	esperado[6] = 0.0;
	esperado[7] = 20488505.595838312;
	esperado[8] = -16734286.967633706;

	double r[POL_DEG];
	int nRoots = roots(p, POL_DEG, r);

	if(verbose)
	{
		printf("Roots:\n");
		printf("N. raices esperado = %d\n", nRootsEsperado);
		printf("N. raices obtenido = %d\n", nRoots);
		printf("Raíz esperada 1 = %.20lf\n", esperado[0]);
		printf("Raíz obtenida 1 = %.20lf\n", r[0]);
		printf("Raíz esperada 2 = %.20lf\n", esperado[1]);
		printf("Raíz obtenida 2 = %.20lf\n", r[1]);
		printf("Raíz esperada 3 = %.20lf\n", esperado[2]);
		printf("Raíz obtenida 3 = %.20lf\n", r[2]);
		printf("Raíz esperada 4 = %.20lf\n", esperado[3]);
		printf("Raíz obtenida 4 = %.20lf\n", r[3]);
		printf("Raíz esperada 5 = %.20lf\n", esperado[4]);
		printf("Raíz obtenida 5 = %.20lf\n", r[4]);
		printf("Raíz esperada 6 = %.20lf\n", esperado[5]);
		printf("Raíz obtenida 6 = %.20lf\n", r[5]);
		printf("Raíz esperada 7 = %.20lf\n", esperado[6]);
		printf("Raíz obtenida 7 = %.20lf\n", r[6]);
		printf("Raíz esperada 8 = %.20lf\n", esperado[7]);
		printf("Raíz obtenida 8 = %.20lf\n", r[7]);
		printf("Raíz esperada 9 = %.20lf\n", esperado[8]);
		printf("Raíz obtenida 9 = %.20lf\n", r[8]);
	}

	assert(nRoots == nRootsEsperado);
	assert(fabs(r[0] - esperado[0]) < 10e-6);
	assert(fabs(r[1] - esperado[1]) < 10e-6);
	assert(fabs(r[2] - esperado[2]) < 10e-6);
	assert(fabs(r[3] - esperado[3]) < 10e-6);
	assert(fabs(r[4] - esperado[4]) < 10e-6);
	assert(fabs(r[5] - esperado[5]) < 10e-6);
	assert(fabs(r[6] - esperado[6]) < 10e-6);
	assert(fabs(r[7] - esperado[7]) < 10e-6);
	assert(fabs(r[8] - esperado[8]) < 10e-6);
	printf("Roots superado!\n");
}

// Aquí va TestUnit

// Aquí va TestSign

void testFix1(bool verbose)
{
	double n = 20.45;
	int esperado = 20;
	int res = fix(n);

	if(verbose)
	{
		printf("Fix1:\n");
		printf("%.20lf ==> %d\n", n, res);
	}

	assert(res == esperado);
	printf("Fix1 superado!\n");
}

void testFix2(bool verbose)
{
	double n = 19.99;
	int esperado = 19;
	int res = fix(n);

	if(verbose)
	{
		printf("Fix2:\n");
		printf("%.20lf ==> %d\n", n, res);
	}

	assert(res == esperado);
	printf("Fix2 superado!\n");
}

void testFix3(bool verbose)
{
	double n = -6.12;
	int esperado = -6;
	int res = fix(n);

	if(verbose)
	{
		printf("Fix3:\n");
		printf("%.20lf ==> %d\n", n, res);
	}

	assert(res == esperado);
	printf("Fix3 superado!\n");
}

int main(){

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

	// Test cross
	printf("Probando cross!\n");
	testCross1(false);
	testCross2(false);
	testCross3(false);
	printf("cross finalizado!\n\n");

	// Test length

	// Test zeros
	printf("Probando zeros!\n");
	testZeros(false);
	printf("zeros finalizado!\n\n");

	// Test prodMatr
	printf("Probando prodMatr!\n");
	testProdMatr1(false);
	//testProdMatr2(true);
	//testProdMatr3(true);
	printf("prodMatr finalizado!\n\n");

	// Test trans
	printf("Probando trans!\n");
	testTrans(false);
	printf("trans finalizado!\n\n");

	// Test det

	// Test roots
	printf("Probando roots!\n");
	testRoots(false);
	printf("roots finalizado!\n\n");

	// Test unit

	// Test sign

	// Test fix
	printf("Probando fix!\n");
	testFix1(false);
	testFix2(false);
	testFix3(false);
	printf("fix finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n\n");
	return 0;
}
