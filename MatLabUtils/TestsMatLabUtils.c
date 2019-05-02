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

void testNorm2()
{
	double vectorB[] = {
		4935037.85913036,
		472703.320202615,
		3999475.70573182,
		};
	double esperadoB = 6.36976082916145030000e+06;

	if(verbose)
	{
		printf("Norm1:\n");
		printf("Esperado B: %.20lf \n", esperadoB);
		printf("Obtenido B: %.20lf \n", norm(vectorB));
	}

	assert(fabs(esperadoB - norm(vectorB)) < 10e-6);
	printf("Norm2 superado!\n");
}

void testNorm3()
{
	double vectorC[] = {
		4909646.95198536,
		687938.936915757,
		3999494.94894739,
		};
	double esperadoC = 6.36976082916144840000e+06;

	if(verbose)
	{
		printf("Norm1:\n");
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

void testLength(bool verbose)
{

}

int main(){

	// Tests norm
	printf("Probando norm!\n");
	testNorm1(0);
	testNorm2(0);
	testNorm3(0);
	printf("norm finalizado!\n\n");

	// Tests dot
	printf("Probando dot!\n");
	testDot1(0);
	testDot2(0);
	testDot3(0);
	printf("dot finalizado!\n\n");

	// Tests length
	printf("Probando length!\n");
	testLength(0);
	printf("length finalizado!\n\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
