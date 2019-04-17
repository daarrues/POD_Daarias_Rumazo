#include "MatLabUtils.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

typedef int bool;
#define true 1
#define false 0

void testDot1(bool verbose)
{
	double v1[] = {6.33886095165153710000e-01,
			-7.73340379778256310000e-01,
			1.15358294325144310000e-02};

	double v2[] = {4.95099033826459660000e+06,
			2.56563116260380720000e+05,
			3.99946534658133240000e+06};

	double expected = 2.98609046501646750000e+06;

	double res = dot(v1,v2,3);
	
	if(verbose){
		printf("Dot1:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}
	
	assert(fabs(res-expected)<10e-12);
	printf("Dot1 superado!\n");
}

void testDot2(bool verbose)
{
	double v1[] = {9.35539825569648990000e-01,
			-1.64830818224054180000e-02,
			-3.52836424971610050000e-01};

	double v2[] = {4.93503785913035740000e+06,
			4.72703320202615400000e+05,
			3.99947570573182080000e+06};

	double expected = 3.19797214063458420000e+06;

	double res = dot(v1,v2,3);
	
	if(verbose){
		printf("Dot2:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}
	
	assert(fabs(res-expected)<10e-12);
	printf("Dot2 superado!\n");
}

void testDot3(bool verbose)
{
	double v1[] = {2.04354223521544110000e+07,
			1.07069944671824550000e+06,
			1.01290549143365120000e+06};

	double v2[] = {1.71964697862374220000e+01,
			-2.65751027611478370000e+03,
			3.73838685080780670000e+03};

	double expected = 1.29265491105025480000e+09;

	double res = dot(v1,v2,3);
	
	if(verbose){
		printf("Dot3:\n");
		printf("Exp.:%.20lf\n",expected);
		printf("Obt.:%.20lf\n",res);
	}
	
	assert(fabs(res-expected)<10e-12);
	printf("Dot3 superado!\n");
}

void testLength(bool 0)
{
	
}

int main(){

	// Tests length
	printf("Probando length!\n");
	testLength(0);
	printf("length finalizado!\n");


	// Tests dot
	printf("Probando dot!\n");
	testDot1(0);
	testDot2(0);
	testDot3(0);
	printf("dot finalizado!\n");

	// Final
	printf("Todos los test superados!\n");
	return 0;
}
