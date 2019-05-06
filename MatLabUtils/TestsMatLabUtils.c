#include "MatLabUtils.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

typedef int bool;
#define true 1
#define false 0

void testNorm1(bool verbose)
{
	double vectorA[] = {
		4950990.33826460,
		256563.116260381,
		3999465.34658133
		};
	double esperadoA = 6.36976082916144930000e+06;

	if(verbose)
	{
		printf("Norm1:\n");
		printf("Esperado A: %.20lf \n", esperadoA);
		printf("Obtenido A: %.20lf \n", norm(vectorA));
	}

	assert(fabs(esperadoA - norm(vectorA)) < 10e-6);
	printf("Norm1 superado!\n");
}

void testNorm2(bool verbose)
{
	double vectorB[] = {
		4935037.85913036,
		472703.320202615,
		3999475.70573182,
		};
	double esperadoB = 6.36976082916145030000e+06;

	if(verbose)
	{
		printf("Norm2:\n");
		printf("Esperado B: %.20lf \n", esperadoB);
		printf("Obtenido B: %.20lf \n", norm(vectorB));
	}

	assert(fabs(esperadoB - norm(vectorB)) < 10e-6);
	printf("Norm2 superado!\n");
}

void testNorm3(bool verbose)
{
	double vectorC[] = {
		4909646.95198536,
		687938.936915757,
		3999494.94894739,
		};
	double esperadoC = 6.36976082916144840000e+06;

	if(verbose)
	{
		printf("Norm3:\n");
		printf("Esperado C: %.20lf \n", esperadoC);
		printf("Obtenido C: %.20lf \n", norm(vectorC));
	}

	assert(fabs(esperadoC - norm(vectorC)) < 10e-6);
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

void testLength(bool verbose)
{

}

int main(){

	// Test norm
	printf("Probando norm!\n");
	testNorm1(0);
	testNorm2(0);
	testNorm3(0);
	printf("norm finalizado!\n\n");

	// Test dot
	printf("Probando dot!\n");
	testDot1(0);
	testDot2(0);
	testDot3(0);
	printf("dot finalizado!\n\n");

	// Tests Cross
	printf("Probando cross!\n");
	testCross1(false);
	testCross2(false);
	testCross3(false);
	printf("cross finalizado!\n\n");

	// Test length
	printf("Probando length!\n");
	testLength(0);
	printf("length finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
