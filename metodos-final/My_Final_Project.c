#include "My_Final_Project.h"

float FuncEval(float x, float y){
    //return (2 * x);
    return 7*ElPow(2, x)*LogNat(2, 16);
}

double Integrador(double x, double p) {
    printf("\tIntegrador activado...");
    double newtonMod = NewtonRapsonMod(x, 5);
    printf("\n\t\tNewton Rapson Mod: %f\n", newtonMod);
    return 0.0;

}

matriz RoungeKoutta(float xi, float yi, float xf) {
    double k1, k2, k3, k4;
    matriz ans = {20, 2}; //Iniciamos matriz con xf + 1 para capturar el ultimo valor de la iteracion
    float h = (xf - xi) / 19;
    ans.mtx[0][0] = xi; //Ingresamos valores iniciales a tabla
    ans.mtx[0][1] = yi;
    int i = 1; //contador iniciando en 1
    while (xi < xf) {
        k1 = FuncEval(xi, yi);
        k2 = FuncEval(xi + h / 2, yi + k1 * h / 2);
        k3 = FuncEval(xi + h / 2, yi + k2 * h /2);
        k4 = FuncEval(xi + h, yi + k3 * h);
        xi = xi + h;
        yi = yi + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
        //printf("\n\t xi: %f and yi %f", xi, yi);
        ans.mtx[i][0] = xi; //ingresamos valores despues de ser modificados a tabla
        ans.mtx[i][1] = yi;
        i++;//se suma a la iteracion
    }
    return(ans);
}

float Correlacion(matriz mat) {
 //   printf("\n\n\tFuncion de correlacion activada...");
    int i;
    float sumatoriaXY, sumatoriaXi, sumatoriaX2, sumatoriaYi, sumatoriaY2, ans;
    sumatoriaX2 = sumatoriaXi = sumatoriaXY = sumatoriaYi = sumatoriaY2 = ans = 0.0;
    float divisor = 0.0, dividendo = 0.0;
    for(i = 0; i < mat.ren; i++) {
        sumatoriaXi += mat.mtx[i][0];
        sumatoriaYi += mat.mtx[i][1];
        sumatoriaX2 += MyPow(mat.mtx[i][0], 2);
        sumatoriaY2 += MyPow(mat.mtx[i][1], 2);
        sumatoriaXY += mat.mtx[i][0] * mat.mtx[i][1];
    }
    divisor = mat.ren * sumatoriaXY - sumatoriaXi * sumatoriaYi;
    dividendo = sqrt(mat.ren * sumatoriaX2 - MyPow(sumatoriaXi, 2)) * sqrt((mat.ren * sumatoriaY2 - MyPow(sumatoriaYi, 2)));
    ans = divisor /dividendo;
 //   printf("\n\t Correlacion es: %f", ans);
    return(ans);
}

void BestFit(matriz t, float x) {
    matriz Array = {7, 1};
    int i = 0;
    int idx = 0;
    float myNumber = 1.0;
    float theNumber = 0.0;
    float linearReg =  Array.mtx[0][0]=  RegLinBestFit(t); //Lineal -> a0 + a1X -
    float cuadraticReg = Array.mtx[1][0] = RegCuadBestFit(t,x);//Cuadratica -> a0 + a1X + a2X^2 x -
    float cubicReg = Array.mtx[2][0] = RegCubicBestFit(t, x);//Cubica -> a0 + a1X + a2X^2 + a3X^3 x -
    float logaritmicReg = Array.mtx[4][0] = RegLnBestFit(t);//Logaritmica -> a0 + a1 * ln(x) -
    float inverseReg = Array.mtx[5][0] = Abs(RegInvBestFit(t));//Inversa -> a0 + a1/x
    float expontialReg = Array.mtx[3][0] = RegExpBestFit(t);//Exponencial -> a0 + e^(A1X) -
    float potenciaReg = Array.mtx[6][0] =RegPowBestFit(t);//Potencia -> a0*X^a1
    float distance = Abs(Array.mtx[0][0] - myNumber);
    for(i; i < Array.ren; i++) {
        float idistance = Abs(Array.mtx[i][0] - myNumber);
        if(idistance < distance){
            idx = i;
            distance = idistance;
        }
    }
    theNumber = Array.mtx[idx][0];
    PrintMtx(Array);
    RegresionChecker(idx, t, x);
}

void RegresionChecker(int idx, matriz t, float x) {
    switch(idx) {
        case 0:
            RegLin(t, x);
            break;
        case 1:
            RegCuad(t, x);
            break;
         case 2:
            RegCubica(t, x);
            break;
        case 3:
            RegExp(t, x);
            break;
        case 4:
            RegLin(t, x);
            break;
        case 5:
            RegInv(t, x);
            break;
        case 6:
            RegPow(t, x);
            break;
    }
}

void Interpolador(float xi, float yi, float h, float xf) {
    matriz tabla = RoungeKoutta(xi, yi, xf);
    //PrintMtx(tabla);
    BestFit(tabla, 1);
    //Correlacion(tabla);
}

float RegLinBestFit(matriz t){
    return Correlacion(t);
}

float RegCuadBestFit(matriz t, float x) {
    matriz regresion = RegCuadMat(t, x);
    matriz ans = t;
    int i = 0;
    for(i; i < t.ren; i++) {
        ans.mtx[i][0] = regresion.mtx[0][0] + regresion.mtx[1][0] * t.mtx[i][0] + regresion.mtx[2][0] * MyPow(t.mtx[i][0] ,2);
    }
    return Correlacion(ans);
}

float RegCubicBestFit(matriz t, float x) {
    matriz regresion = RegCubicMat(t, x);
    matriz ans = t;
    int i = 0;
    for(i; i < t.ren; i++) {
        ans.mtx[i][0] = regresion.mtx[0][0] + regresion.mtx[1][0] * t.mtx[i][0] + regresion.mtx[2][0] * MyPow(t.mtx[i][0] ,2) + regresion.mtx[3][0] * MyPow(t.mtx[i][0], 3);
    }
    return Correlacion(ans);
}

float RegExpBestFit(matriz t) {
    matriz ans = t;
    printf("\n\tVerificacion de regresion Exponensial favor de esperar un momento...\n\t");
    int i = 0;
    for(i; i < t.ren; i++){
        ans.mtx[i][1] = LogNat(t.mtx[i][1], 5);
    }
    return Correlacion(ans);
}

float RegLnBestFit(matriz t) {
    matriz ans = t;
    int i = 0;
    for(i; i < t.ren; i++){
        ans.mtx[i][0] = LogNat(t.mtx[i][0], 5);
    }
    return Correlacion(ans);
}

float RegInvBestFit(matriz t) {
    matriz ans = t;
    int i = 0;
    for(i; i < t.ren; i++){
        ans.mtx[i][0] = 1/(t.mtx[i][0]);
    }
    return Correlacion(ans);
}

float RegPowBestFit(matriz t) {
    matriz ans = t;
    int i = 0;
    printf("\n\tVerificacion de regresion de Potencia favor de esperar un momento...\n\t");
    for(i; i < t.ren; i++){
        ans.mtx[i][0] = LogNat(t.mtx[i][0], 5);
        ans.mtx[i][1] = LogNat(t.mtx[i][1], 5);
    }
    return Correlacion(ans);
}

matriz RegCuadMat(matriz t, float x) {
    float ans = 0.0;
    matriz a;
    a = Regresion(t, 2);
    return (a);
}

matriz RegCubicMat(matriz t, float x) {
    float ans = 0.0;
    matriz a;
    a = Regresion(t, 3);
    ans = a.mtx[0][0] + a.mtx[1][0] * x + a.mtx[2][0], 2 * MyPow(x, 2) + a.mtx[3][0] * MyPow(x, 3);
    return (a);
  }


