#include "MatLabUtils.h"
#include <stdio.h>
#define size(v)  (sizeof(v) / sizeof((v)[0]))

double norm(double* v)
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double dot(double *v1, double *v2, int dim)
{
	double res = 0;
	for(int i = 0; i < dim; i++)
	{
		res += v1[i]*v2[i];
	}
	return res;
}

double length(double *v)
{
	return size(v);
}

double[][] prodMatr(double[][] m1, double[][] m2)
{
	for(int j = 0; j < 3; j++){

	}
}
