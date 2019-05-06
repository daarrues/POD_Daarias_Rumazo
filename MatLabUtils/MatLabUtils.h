#ifndef MATLABUTILS_H
#define MATLABUTILS_H

double norm(double[3]);
double dot(double[3], double[3]);
// cross (producto vectorial)
double length(double*);
// zeros 		(inicializar a 0)
void prodMatr(double[3][3], double[3][3], double[3][3]);
// ‘ 		(traspuesta)
double det(double[3][3]);
// roots		(raíces de polinomio)
void unit(double[3]);
// isreal		(pertenencia a los reales) // <-- QUIZAS NO HAGA FALTA POR ROOTS?
int sign(double);
// fix

#endif
