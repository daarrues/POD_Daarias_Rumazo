#ifndef MATLABUTILS_H
#define MATLABUTILS_H

double norm(double[3]);
double dot(double[3], double[3]);
void cross(double[3], double[3], double[3]);
double length(double*);
void zeros(double[][], int rows, int cols);
void prodMatr(double[3][3], double[3][3], double[3][3]);
<<<<<<< HEAD
// ‘ 		(traspuesta)
double det(double[3][3]);
=======
void trans(double[3][3], double[3][3]);
// det 		(determinante)
>>>>>>> a52bb767695ab5d48a0007ca782ca96ccdbda66d
// roots		(raíces de polinomio)
void unit(double[3]);
// isreal		(pertenencia a los reales) // <-- QUIZAS NO HAGA FALTA POR ROOTS?
int sign(double);
// fix

#endif
