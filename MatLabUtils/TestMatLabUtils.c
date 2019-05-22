//------------------------------------------------------------------------------
//                              TestMatLabUtils
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file TestMatLabUtils.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/04/17
 *
 * Este fichero contiene los test unitarios
 * para las funciones de alto nivel de MatLab
 * necesarias para este proyecto.
 */
//------------------------------------------------------------------------------
#include "MatLabUtils.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

typedef int bool;
#define true 1
#define false 0

//------------------------------------------------------------------------------
//  void testNorm1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función norm
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testNorm2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función norm
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testNorm3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función norm
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------
//  void testDot1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función dot
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testDot2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función dot
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testDot3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función dot
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testCross1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función cross
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testCross2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función cross
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testCross3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función cross
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testZeros(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación de la función zeros
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testZeros(bool verbose)
{
	double v[3];
	zeros(v, 3);

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

//------------------------------------------------------------------------------
//  void testProdMatr1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función prodMatr
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testProdMatr2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función prodMatr
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testProdMatr2(bool verbose)
{
	// Matriz 1
	double m1[3][3] = {
		{
			0.999999997895568,
			-5.95228911876096e-05,
			-2.58048231218136e-05
		},
		{
			5.95223155942333e-05,
			0.999999997979772,
			-2.23058450876347e-05
		},
		{
			2.5806150778072e-05,
			2.23043090778784e-05,
			0.99999999941828
		}
	};

	// Matriz 2
	double m2[3][3] = {
		{
			0.999997373791705,
			-0.00210195244551542,
			-0.000913348569948463
		},
		{
			0.00210195244550126,
			0.999997790895058,
			-9.59924398978515e-07
		},
		{
			0.000913348569981061,
			-9.59893382181152e-07,
			0.999999582896647
		}
	};

	// Matriz esperada
	double mE[3][3] = {
		{
			0.999997223004194,
			-0.00216147517601742,
			-0.000939153323247443
		},
		{
			0.00216145422751936,
			0.999997663783168,
			-2.33201248026592e-05
		},
		{
			0.000939201535052504,
			2.12901231219589e-05,
			0.999999558723506
		}
	};

	// Matriz resultado
	double mR[3][3];

	prodMatr(m1, m2, mR);

	if(verbose) printf("ProdMatr2:\n");
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
	printf("ProdMatr2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testProdMatr3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función prodMatr
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testProdMatr3(bool verbose)
{
				// Matriz 1
				double m1[3][3] = {
					{
						0.999999997895152,
						-5.95287721488282e-05,
						-2.58073726784525e-05
					},
					{
						5.95281964998585e-05,
						0.999999997979423,
						-2.23057957973194e-05
					},
					{
						2.58087004629423e-05,
						2.2304259484005e-05,
						0.999999999418215
					}
				};

				// Matriz 2
				double m2[3][3] = {
					{
						0.99999737378108,
						-0.0021019566973475,
						-0.000913350417381566
					},
					{
						0.00210195669733333,
						0.99999779088612,
						-9.59928282397388e-07
					},
					{
						0.000913350417414164,
						-9.59897265411807e-07,
						0.99999958289496
					}
				};

				// Matriz esperada
				double mE[3][3] = {
					{
						0.999997222978162,
						-0.00216148530879376,
						-0.000939157720229817
					},
					{
						0.00216146436024405,
						0.999997663761267,
						-2.33200848770364e-05
					},
					{
						0.000939205932154037,
						2.12900641757145e-05,
						0.999999558719378
					}
				};

				// Matriz resultado
				double mR[3][3];

				prodMatr(m1, m2, mR);

				if(verbose) printf("ProdMatr3:\n");
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
				printf("ProdMatr3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testTrans(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación de la función trans
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testDet(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación de la función det
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testDet(bool verbose)
{
	double m[3][3] = {
		{
			0.953375083066611,
			0.976292565505083,
			0.991769286600702
		},
		{
			0.21677824571211,
			0.163060061341653,
			0.106138375728988
		},
		{
			-0.209959860863242,
			-0.142352530482458,
			-0.0716123407880656
		}
	};

	double expected = 2.075911e-05;

	double res = det(m);

	if(verbose)
	{
		printf("Det:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}

	assert(fabs(res-expected)<10e-12);
	printf("Det superado!\n");
}

//------------------------------------------------------------------------------
//  void testRoots(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación de la función roots
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testUnit1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función unit
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testUnit1(bool verbose)
{
	double v[3] = {
		20408482.6529284,
		1869905.83442928,
		-114536.67089787
	};

	double vE[3] = {
		0.995813218529303,
		0.0912403425083912,
		-0.0055887119501299
	};

	double vR[3];

	unit(v, vR);

	if(verbose) printf("Unit1:\n");
	for(int i = 0; i < 3; i++)
	{
			if(verbose)
			{
				printf("Exp.:%.20lf\n", vE[i]);
				printf("Obt.:%.20lf\n", vR[i]);
			}
			assert(fabs(vR[i]-vE[i])<10e-12);
	}
	printf("Unit1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testUnit2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función unit
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testUnit2(bool verbose)
{
	double v[3] = {
		2003323683191.96,
		-22852239115271.2,
		-16251616553379
	};

	double vE[3] = {
		0.0712591311936418,
		-0.812864500552877,
		-0.578077365032094
	};

	double vR[3];

	unit(v, vR);

	if(verbose) printf("Unit2:\n");
	for(int i = 0; i < 3; i++)
	{
			if(verbose)
			{
				printf("Exp.:%.20lf\n", vE[i]);
				printf("Obt.:%.20lf\n", vR[i]);
			}
			assert(fabs(vR[i]-vE[i])<10e-12);
	}
	printf("Unit2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testUnit3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función unit
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testUnit3(bool verbose)
{
	double v[3] = {
		20370508.3427397,
		1861271.24297514,
		-106173.66558631
	};

	double vE[3] = {
		0.995838223643539,
		0.0909906133483308,
		-0.00519043475775007
	};

	double vR[3];

	unit(v, vR);

	if(verbose) printf("Unit3:\n");
	for(int i = 0; i < 3; i++)
	{
			if(verbose)
			{
				printf("Exp.:%.20lf\n", vE[i]);
				printf("Obt.:%.20lf\n", vR[i]);
			}
			assert(fabs(vR[i]-vE[i])<10e-12);
	}
	printf("Unit3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testSign1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función sign
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testSign1(bool verbose)
{
	double n = 34.2348293234;
	int esperado = 1;
	int res = sign(n);

	if(verbose)
	{
		printf("Sign1:\n");
		printf("%.20lf ==> %d\n", n, res);
	}

	assert(res == esperado);
	printf("Sign1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testSign2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función sign
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testSign2(bool verbose)
{
	double n = -12.0905823;
	int esperado = -1;
	int res = sign(n);

	if(verbose)
	{
		printf("Sign2:\n");
		printf("%.20lf ==> %d\n", n, res);
	}

	assert(res == esperado);
	printf("Sign2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testSign3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función sign
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testSign3(bool verbose)
{
	double n = 0.0;
	int esperado = 0;
	int res = sign(n);

	if(verbose)
	{
		printf("Sign3:\n");
		printf("%.20lf ==> %d\n", n, res);
	}

	assert(res == esperado);
	printf("Sign3 superado!\n");
}

//------------------------------------------------------------------------------
//  void testFix1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función fix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testFix2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función fix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testFix3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función fix
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void testAll1(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 1 de la función all
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testAll1(bool verbose)
{
	double v[] = {
		4950990.33826460,
		256563.116260381,
		3999465.34658133
		};
	int esperado = 0;
	int obtenido = all(v);

	if(verbose)
	{
		printf("All1:\n");
		printf("Esperado: %d \n", esperado);
		printf("Obtenido: %d \n", obtenido);
	}

	assert(esperado == obtenido);
	printf("All1 superado!\n");
}

//------------------------------------------------------------------------------
//  void testAll2(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 2 de la función all
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testAll2(bool verbose)
{
	double v[] = {
		0.0,
		0.000001,
		0.0
		};
	int esperado = 0;
	int obtenido = all(v);

	if(verbose)
	{
		printf("All2:\n");
		printf("Esperado: %d \n", esperado);
		printf("Obtenido: %d \n", obtenido);
	}

	assert(esperado == obtenido);
	printf("All2 superado!\n");
}

//------------------------------------------------------------------------------
//  void testAll3(bool verbose)
//------------------------------------------------------------------------------
/**
 * Comprobación 3 de la función all
 *
 * @param <verbose> booleano que indica si debe mostrar lo esperado y obtenido.
 */
//------------------------------------------------------------------------------
void testAll3(bool verbose)
{
	double v[] = {
		0.0,
		0.0,
		0.0
		};
	int esperado = 1;
	int obtenido = all(v);

	if(verbose)
	{
		printf("All3:\n");
		printf("Esperado: %d \n", esperado);
		printf("Obtenido: %d \n", obtenido);
	}

	assert(esperado == obtenido);
	printf("All3 superado!\n");
}

//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
 * Función principal, realiza todas las comprobaciones.
 */
//------------------------------------------------------------------------------
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

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
