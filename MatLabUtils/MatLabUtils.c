//------------------------------------------------------------------------------
//                              MatLabUtils
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file MatLabUtils.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/04/17
 *
 * Este fichero contiene las implementaciones
 * para las funciones de alto nivel de MatLab
 * necesarias para este proyecto.
 */
//------------------------------------------------------------------------------
#include "MatLabUtils.h"
#include "rpoly/rpoly.h"
#include <stdio.h>
#include <math.h>

#define POL_DEG 15

//------------------------------------------------------------------------------
//  double norm(double v[3])
//------------------------------------------------------------------------------
/**
 * Calcula la norma de un vector de 3 componentes reales.
 *
 * @param <v> vector de 3 componentes reales.
 *
 * @return Norma del vector.
 */
//------------------------------------------------------------------------------
double norm(double v[3])
{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

//------------------------------------------------------------------------------
//  double dot(double v1[3], double v2[3])
//------------------------------------------------------------------------------
/**
 * Calcula el producto escalar de dos vectores de 3
 * componentes reales
 *
 * @param <v1> vector de 3 componentes reales.
 * @param <v2> vector de 3 componentes reales.
 *
 * @return Producto escalar de v1 y v2.
 */
//------------------------------------------------------------------------------
double dot(double v1[3], double v2[3])
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

//------------------------------------------------------------------------------
//  void cross(double v1[3], double v2[3], double vResult[3])
//------------------------------------------------------------------------------
/**
 * Calcula el producto vectorial de dos vectores de 3 componentes reales
 *
 * @param <v1>      primer vector de entrada.
 * @param <v2>      segundo vector de entrada.
 * @param <vResult> vector de salida para almacenar el producto escalar.
 */
//------------------------------------------------------------------------------
void cross(double v1[3], double v2[3], double vResult[3])
{
	vResult[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vResult[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vResult[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

//------------------------------------------------------------------------------
//  void zeros(double v[], int n)
//------------------------------------------------------------------------------
/**
 * Inicializa a 0 las componentes de un vector real.
 *
 * @param <v> vector a inicializar a 0.
 * @param <n> número de componentes del vector.
 */
//------------------------------------------------------------------------------
void zeros(double v[], int n)
{
	for(int i = 0; i < n; i++)
	{
		v[i] = 0.0;
	}
}

//------------------------------------------------------------------------------
//  double prodMatr(double m1[3][3], double m2[3][3], double mResult[3][3])
//------------------------------------------------------------------------------
/**
 * Calcula el producto de dos matrices de dimensión 3x3
 * de componentes reales
 *
 * @param <m1>      matriz 3x3 de componentes reales.
 * @param <m2>      matriz 3x3 de componentes reales.
 * @param <mResult> matriz de salida para almacenar el producto de m1 por m2.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void trans(double m[3][3], double mResult[3][3])
//------------------------------------------------------------------------------
/**
 * Calcula la matriz traspuesta de una matriz cuadrada de orden 3.
 *
 * @param <m> matriz de entrada.
 * @param <t> matriz de salida para almacenar la traspuesta de m.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  double det(double m[3][3])
//------------------------------------------------------------------------------
/**
 * Calcula el determinante de una matriz de dimensión 3x3
 * de componentes reales
 *
 * @param <m> matriz 3x3 de componentes reales.
 *
 * @return Determinante de m.
 */
//------------------------------------------------------------------------------
double det(double m[3][3])
{
	return m[0][0] * m[1][1] * m[2][2] +
				 m[0][1] * m[1][2] * m[2][0] +
				 m[1][0] * m[0][2] * m[2][1] -
				 m[0][2] * m[1][1] * m[2][0] -
				 m[0][1] * m[1][0] * m[2][2] -
				 m[0][0] * m[1][2] * m[2][1];
}

//------------------------------------------------------------------------------
//  int roots(double p[], int degree, double r[])
//------------------------------------------------------------------------------
/**
 * Encuentra las raíces reales de un polinomio con coeficientes reales.
 *
 * @param <p>      polinomio con coeficientes reales.
 * @param <degree> grado del polinomio.
 * @param <r>      vector de salida para almacenar las raíces reales de p.
 *
 * @return Número de raíces reales de p incluyendo su multiplicidad.
 */
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//  void unit(double v[3], double vR[3])
//------------------------------------------------------------------------------
/**
 * Calcula el vector unitario a partir de
 * uno de 3 componentes reales
 *
 * @param <v> vector de 3 componentes reales.
 * @param <vR> vector de salida para almacenar el vector unitario respecto a v.
 */
//------------------------------------------------------------------------------
void unit(double v[3], double vR[3])
{
	double n = norm(v);
	vR[0] = v[0]/n;
	vR[1] = v[1]/n;
	vR[2] = v[2]/n;
}

//------------------------------------------------------------------------------
//  int sign(double n)
//------------------------------------------------------------------------------
/**
 * Calcula el signo de un número real
 *
 * @param <n> número real.
 *
 * @return -1 si es negativo, 0 si es cero, 1 si es positivo.
 */
//------------------------------------------------------------------------------
int sign(double n)
{
	return (n == 0) ? 0 : (n > 0) ? 1 : -1;
}

//------------------------------------------------------------------------------
//  int fix(double n)
//------------------------------------------------------------------------------
/**
 * Redondea hacia el 0 (trunca) un número real.
 * Devuelve el mayor entero por debajo para valores positivos (similar a floor),
 * y el menor entero por encima para valores negativos (similar a ceil).
 *
 * @param <n> número real.
 *
 * @return n truncado a los enteros.
 */
//------------------------------------------------------------------------------
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
