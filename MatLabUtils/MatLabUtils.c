#include "MatLabUtils.h"
#include <stdio.h>
#include <math.h>
#define size(v)  (sizeof(v) / sizeof((v)[0]))

double norm(double v[3])
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double dot(double v1[3], double v2[3])
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double length(double *v)
{
	return size(v);
}

void prodMatr(double m1[3][3], double m2[3][3], double mResult[3][3])
{
	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			mResult[i][j] = dot(m1[i],m2[j]);
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

void unit(double v[3])
{
	double norm = norm(v);
	v[0] = v[0]/norm;
	v[1] = v[1]/norm;
	v[2] = v[2]/norm;
}

int sign(double n)
{
	return (n == 0) ? 0 : (n > 0) ? 1 : -1;
}
