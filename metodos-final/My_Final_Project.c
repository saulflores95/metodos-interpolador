#include "My_Final_Project.h"
#include "My_Intergral_dif.h"

float FuncEval(float x){
    return (2 * x);
}

matriz RoungeKoutta(float xi, float yi, float h, float xf) {
    double k1, k2, k3, k4;
    matriz ans = {xf + 1, 2}; //Iniciamos matriz con xf + 1 para capturar el ultimo valor de la iteracion
    ans.mtx[0][0] = xi; //Ingresamos valores iniciales a tabla
    ans.mtx[0][1] = yi;
    int i = 1; //contador iniciando en 1
    while (xi < xf) {
        k1 = Fn(xi, yi);
        k2 = Fn(xi + h / 2, yi + k1 * h / 2);
        k3 = Fn(xi + h / 2, yi + k2 * h /2);
        k4 = Fn(xi + h, yi + k3 * h);
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
    printf("\n\tFuncion de correlacion activada...");
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
    PrintTable(t);
    //Lineal -> a0 + a1X
    printf("\n\t x = %f", x);
    RegLin(t, x);
    //printf("\n\t La evaluacion de Y(%f) = %f", x, ans);
    //Cuadratica -> a0 + a1X + a2X^2
    //Cubica -> a0 + a1X + a2X^2 + a3X^3
    //Exponencial -> a0 + e^(A1X)
    //Logaritmica -> a0 + a1 * ln(x)
    //Inversa -> a0 + a1/x
    //Potencia -> a0*X^a1

}

void Interpolador(float xi, float yi, float h, float xf) {
    matriz tabla = RoungeKoutta(xi, yi, h, xf);
    PrintMtx(tabla);
    //BestFit(tabla, 4.0);
    Correlacion(tabla);
}
