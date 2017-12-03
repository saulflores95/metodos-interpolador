#include "My_Intergral_dif.h"
#include <math.h>
//Funcion a intergar
double FI(float x){
    double ans = 0.0;
    //ans = 400 * MyPow(x, 5) - 900*MyPow(x,4) + 675 * MyPow(x, 3) - 200 * MyPow(x,2) + 25 * x + 0.2;
    ans = exp(-x * x / 2.0)/ sqrt(PI);
    return(ans);
}
//Metodo del trapecio multiple
double TrapecioM(float a, float b, int n) {
    double ans = 0.0, h, sum = 0.0;
    int i;
    h = (b - a)/n;
    for(i = 1; i < n; i++)
        sum += FI(a + i * h);
    ans = (FI(a) + 2.0 * sum + FI(b)) * h / 2.0;
    return(ans);
}
//Regla de Simpson 1/3 aplicacion multiple
double Simpson3M(float a, float b, int n){
    double ans = 0.0, sp = 0.0, si = 0.0, h; //sp = suma par y si = suma impar
    int i;
    if(Residuo(n, 2.0) != 0.0)
        n++;
    h = (b - a) / n;
    for(i = 1; i < n - 1; i += 2)
        si += FI(a + i * h);
    for(i = 2; i < n - 2; i+= 2)
        sp += FI(a + i * h);
    ans = (FI(a) + 4.0 * si + 2.0 * sp + FI(b) * h / 3.0);
    return(ans);
}
//Metodo rounge kutta de cuarto orden
double Fn(float x, float y) {
    double ans = 0.0;
    ans = 4 * exp(0.8 * x) - 0.5 * y;
    return(ans);
}

double R4kc(float xi, float yi, float h, float xf) {
    double k1, k2, k3, k4;
    while (xi < xf) {
        k1 = Fn(xi, yi);
        k2 = Fn(xi + h / 2, yi + k1 * h / 2);
        k3 = Fn(xi + h / 2, yi + k2 * h /2);
        k4 = Fn(xi + h, yi + k3 * h);
        xi = xi + h;
        yi = yi + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
        printf("\n\t xi: %f and yi %f", xi, yi);
    }
    return(yi);
}
