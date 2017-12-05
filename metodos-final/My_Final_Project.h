#ifndef MY_FINAL_PROJECT_H_INCLUDED
#define MY_FINAL_PROJECT_H_INCLUDED
#include "My_Intergral_dif.h"
#include "My_Reg_Inter.h"

matriz RoungeKoutta(float xi, float yi, float xf);

double Integrador(double, double);

float FuncEval(float, float);

float Correlacion(matriz);

void RegresionChecker(int, matriz, float);

void BestFit(matriz, float);

void Interpolador(float xi, float yi, float h, float xf);

float RegLinBestFit(matriz);

float RegCuadBestFit(matriz, float);

float RegCubicBestFit(matriz, float);

float RegExpBestFit(matriz t);

float RegLnBestFit(matriz t);

float RegInvBestFit(matriz);

float RegPowBestFit(matriz t);

matriz RegCuadMat(matriz, float);

matriz RegCubicMat(matriz, float);
#endif // MY_FINAL_PROJECT_H_INCLUDED
