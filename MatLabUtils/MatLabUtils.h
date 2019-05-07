#ifndef MATLABUTILS_H
#define MATLABUTILS_H

double norm(double[3]);
double dot(double[3], double[3]);
void cross(double[3], double[3], double[3]);
double length(double*);
void zeros(double[3]);
void prodMatr(double[3][3], double[3][3], double[3][3]);
void trans(double[3][3], double[3][3]);
double det(double[3][3]);
int roots(double[], int, double[]);
void unit(double[3]);
// isreal		(pertenencia a los reales) // <-- QUIZAS NO HAGA FALTA POR ROOTS?
int sign(double);
int fix(double);

#endif
