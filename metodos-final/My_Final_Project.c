#include "My_Final_Project.h"
#include "My_Intergral_Dif.h"

float FuncEval(float x, float y){
    return (2 * x);
}

double Integrador(double P, double x) {
    printf("Integrador activado...");
    double newton = NewtonRapson(0.5, 5);
    //printf("\nNewton: %f", newton);
    printf("\nSimpson %f", Simpson3M(-10, 5, 500));
    return (Simpson3M(-10, 5, 500) - P);

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
    printf("\n\n\tFuncion de correlacion activada...");
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
    printf("\n\t Correlacion es: %f", ans);
    return(ans);
}

void BestFit(matriz t, float x) {
    /*
    float linearReg = RegLinBestFit(t); //Lineal -> a0 + a1X -
    float cuadraticReg = RegCuadBestFit(t,x);//Cuadratica -> a0 + a1X + a2X^2 x -
    float cubicReg = RegCubicBestFit(t, x);//Cubica -> a0 + a1X + a2X^2 + a3X^3 x -
    float expontialReg = RegExpBestFit(t);//Exponencial -> a0 + e^(A1X) -
    float logaritmicReg = RegLnBestFit(t);//Logaritmica -> a0 + a1 * ln(x) -
    float inverseReg = RegInvBestFit(t);//Inversa -> a0 + a1/x
    float potenciaReg = RegPowBestFit(t);//Potencia -> a0*X^a1
    */
}

void Interpolador(float xi, float yi, float h, float xf) {
    matriz tabla = RoungeKoutta(xi, yi, xf);
    PrintMtx(tabla);
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
    printf("\n\tLa regresion cuadratica es:");

    if (a.mtx[1][0] >= 0.0 && a.mtx[2][0] >= 0.0)
      printf("\n\t y = %f + %f x + %f x^2", a.mtx[0][0], a.mtx[1][0], a.mtx[2][0]);

    else if (a.mtx[1][0] < 0.0 && a.mtx[2][0] >= 0.0)
      printf("\n\t y = %f - %f x + %f x^2", a.mtx[0][0], -a.mtx[1][0], a.mtx[2][0]);

    else if (a.mtx[1][0] >= 0.0 && a.mtx[2][0] < 0.0)
      printf("\n\t y = %f + %f x - %f x^2", a.mtx[0][0], a.mtx[1][0], -a.mtx[2][0]);

    else if (a.mtx[1][0] < 0.0 && a.mtx[2][0] < 0.0)
      printf("\n\t y = %f - %f x - %f x^2", a.mtx[0][0], -a.mtx[1][0], -a.mtx[2][0]);

    ans = a.mtx[0][0] + a.mtx[1][0] * x + a.mtx[2][0], 2 * MyPow(x, 2);
    printf("\n\t La evaluacion de Y(%f) = %f", x, ans);

    return (a);
}

matriz RegCubicMat(matriz t, float x) {
    float ans = 0.0;
    matriz a;
    a = Regresion(t, 3);
    printf("\n\tLa regresion cubica es:");
    if (a.mtx[1][0] >= 0.0 && a.mtx[2][0] >= 0.0 && a.mtx[3][0] >= 0.0) // + + +
      printf("\n\t y = %f + %f x + %f x^2 + %f x^3", a.mtx[0][0], a.mtx[1][0], a.mtx[2][0], a.mtx[3][0]);

    else if (a.mtx[1][0] >= 0.0 && a.mtx[2][0] >= 0.0 && a.mtx[3][0] < 0.0) // + + -
      printf("\n\t y = %f + %f x + %f x^2 - %f x^3", a.mtx[0][0], a.mtx[1][0], a.mtx[2][0], -a.mtx[3][0]);

    else if (a.mtx[1][0] >= 0.0 && a.mtx[2][0] < 0.0 && a.mtx[3][0] >= 0.0) // + - +
      printf("\n\t y = %f + %f x - %f x^2 + %f x^3", a.mtx[0][0], a.mtx[1][0], -a.mtx[2][0], a.mtx[3][0]);

    else if (a.mtx[1][0] >= 0.0 && a.mtx[2][0] < 0.0 && a.mtx[3][0] < 0.0) // + - -
      printf("\n\t y = %f + %f x - %f x^2 - %f x^3", a.mtx[0][0], a.mtx[1][0], -a.mtx[2][0], -a.mtx[3][0]);

    else if (a.mtx[1][0] < 0.0 && a.mtx[2][0] >= 0.0 && a.mtx[3][0] >= 0.0) // - + +
      printf("\n\t y = %f - %f x + %f x^2 + %f x^3", a.mtx[0][0], -a.mtx[1][0], a.mtx[2][0], a.mtx[3][0]);

    else if (a.mtx[1][0] < 0.0 && a.mtx[2][0] >= 0.0 && a.mtx[3][0] < 0.0) // - + -
      printf("\n\t y = %f - %f x + %f x^2 - %f x^3", a.mtx[0][0], -a.mtx[1][0], a.mtx[2][0], -a.mtx[3][0]);

    else if (a.mtx[1][0] < 0.0 && a.mtx[2][0] < 0.0 && a.mtx[3][0] >= 0.0) // - - +
      printf("\n\t y = %f - %f x - %f x^2 + %f x^3", a.mtx[0][0], -a.mtx[1][0], -a.mtx[2][0], a.mtx[3][0]);

    else if (a.mtx[1][0] < 0.0 && a.mtx[2][0] < 0.0 && a.mtx[3][0] < 0.0) // - - -
      printf("\n\t y = %f - %f x - %f x^2 - %f x^3", a.mtx[0][0], -a.mtx[1][0], -a.mtx[2][0], -a.mtx[3][0]);

    ans = a.mtx[0][0] + a.mtx[1][0] * x + a.mtx[2][0], 2 * MyPow(x, 2) + a.mtx[3][0] * MyPow(x, 3);
    printf("\n\t La evaluacion de Y(%f) = %f", x, ans);
    return (a);
  }


