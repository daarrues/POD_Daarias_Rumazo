#ifndef MATLABUTILS_H
#define MATLABUTILS_H

double norm(double[3]);
double dot(double[3], double[3]);
void cross(double[3], double[3], double[3]);
double length(double*);
void zeros(double[][], int rows, int cols);
void prodMatr(double[3][3], double[3][3], double[3][3]);
void trans(double[3][3], double[3][3]);
// det 		(determinante)
// roots		(raíces de polinomio)
// unit		(normalizar un vector)
// isreal		(pertenencia a los reales)
// sign		(signo)
// fix

#endif
