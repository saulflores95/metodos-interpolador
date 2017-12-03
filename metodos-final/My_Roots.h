#ifndef MY_ROOTS_H_INCLUDED
#define MY_ROOTS_H_INCLUDED
#include "My_Math.h"
//Funcion a estimar raiz de paracaidista
double FParacaidista(double);

double FExamen(double);
//funcion a estimar raiz
double F(double);
//Derivada a la funcion a estimar raiz
double FD(double);
//funcion que mapea una exprecion
void Tabla(float, float, int);
//funcion que mapea biseccion
void TablaTarea(double, double, int);
//metodo de biseccion
double Biseccion(double, double, int);
//metodo de newton rapson
double NewtonRapson(double, int);

double FExamen(double x);

double FalsaPosicion(float, float, int);

#endif // MY_ROOTS_H_INCLUDED
