#ifndef MY_FINAL_PROJECT_H_INCLUDED
#define MY_FINAL_PROJECT_H_INCLUDED
#include "My_Intergral_dif.h"
#include "My_Reg_Inter.h"

matriz RoungeKoutta(float xi, float yi, float h, float xf);

float FuncEval(float);

float Correlacion(matriz);

void BestFit(matriz, float);

void Interpolador(float xi, float yi, float h, float xf);


#endif // MY_FINAL_PROJECT_H_INCLUDED
