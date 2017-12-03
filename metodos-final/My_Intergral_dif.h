#ifndef MY_INTERGRAL_DIF_H_INCLUDED
#define MY_INTERGRAL_DIF_H_INCLUDED
#include "My_Math.h"
#include "My_Matrix.h"
double FI(float);

double Fn(float, float);

double Trapecio(float, float, int);

//Regla de Simpson 1/3 aplicacion multiple
double Simpson3M(float, float, int);
//Derivacion hacia adelante
//Derivacion hacia atras
//Implementar formulas
double R4kc(float, float, float, float);

#endif // MY_INTERGRAL_DIF_H_INCLUDED
