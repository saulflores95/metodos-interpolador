#ifndef MY_POLI_ROOTS_H_INCLUDED
#define MY_POLI_ROOTS_H_INCLUDED
#include "My_Math.h"

#define GMAX 10 //grado maximo

typedef struct {
  int grado;
  float pol [GMAX];
} polinomio;

//captura un polinomio
polinomio ScanfPol(void);
//imprimer un polinomio
void PrintPol(polinomio p);
//funcion que evalua un polinomio
float PolEval(polinomio, float);
//funcion que deriva un polinomio
polinomio PolDeriv(polinomio);
//Newton-Rapson para polinomios
double NewtonRPol(polinomio, float, int);

polinomio DivisionSint(polinomio, float);

void BirgeVieta(polinomio);


#endif // MY_POLI_ROOTS_H_INCLUDED
