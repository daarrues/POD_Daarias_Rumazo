#include "MatLabUtils.h"
#include "rpoly.h"
#include <stdio.h>
#include <math.h>
#define size(v)  (sizeof(v) / sizeof((v)[0]))
#define POL_DEG 15

double norm(double v[3])
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double dot(double v1[3], double v2[3])
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void cross(double v1[3], double v2[3], double vResult[3])
{
	vResult[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vResult[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vResult[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double length(double *v)
{
	return size(v);
}

void zeros(double v[], int n)
{
	for(int i = 0; i < n; i++)
	{
		v[i] = 0.0;
	}
}

void prodMatr(double m1[3][3], double m2[3][3], double mResult[3][3])
{
	double m2T[3][3];
	trans(m2, m2T);

	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			mResult[i][j] = dot(m1[i], m2T[j]);
		}
	}
}

void trans(double m[3][3], double mResult[3][3])
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			mResult[i][j] = m[j][i];
		}
	}
}

double det(double m[3][3])
{
	return m[0][0] * m[1][1] * m[2][2] +
				 m[0][1] * m[1][2] * m[2][0] +
				 m[1][0] * m[0][2] * m[2][1] -
				 m[0][2] * m[1][1] * m[2][0] -
				 m[0][1] * m[1][0] * m[2][2] -
				 m[0][0] * m[1][2] * m[2][1];
}

int roots(double p[], int degree, double r[])
{
	double zerosR[POL_DEG];
	double zerosI[POL_DEG];
	int n = real_poly_roots(p, degree, zerosR, zerosI);
	int nRoots = 0;
	for(int i = 0; i < degree; i++)
	{
		if(fabs(zerosI[i]) < 10e-12)
		{
			r[nRoots] = zerosR[i];
			nRoots++;
		}
	}
	return nRoots;
}

void unit(double v[3], double vR[3])
{
	double n = norm(v);
	vR[0] = v[0]/n;
	vR[1] = v[1]/n;
	vR[2] = v[2]/n;
}

int sign(double n)
{
	return (n == 0) ? 0 : (n > 0) ? 1 : -1;
}

int fix(double n)
{
	int res;
	if(n >= 0)
	{
		res = floor(n);
	}
	else
	{
		res = ceil(n);
	}
	return res;
}
