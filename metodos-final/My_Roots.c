#include "My_Roots.h"
double GlobalVariable;

//Funcion a estimar raiz de paracaidista
double FParacaidista(double x) {
  float ans, m, g = 9.81, t = 10.0, v = 40.0, c = 15.0;
  m = x;
  ans = m * g * (1.0 - Expo(-c * t / m, 29)) / c - v;
  return (ans);
}
//Funcion a estimar raiz
double Fold(float x) {
  float ans, p = -3000.0, f = 0.0, i = 12.0 / 12.0, n = 5 * 12, b = 150000;
  b = x;
  ans = p * 100.0 * (1.0 - MyPow(1.0 + i / 100, -n)) / i
  + f * MyPow(1.0 + 1 / 100.0, -n) + b;
  return (ans);
}

double F(double x){
  printf("\nIniciado F");
  double ans;
  //ans = LogNat(x, 4) + MyPow(x, 2) - 4;
  //ans = cos(x) - x;
  //ans = (1/sqrt(2 * PI)) * exp(-MyPow(x, 2) / 2);
  ans = Simpson3M(-5, x, 500)- GlobalVariable;
  return (ans);
}
//Derivada a la funcion a estimar raiz
double FD(double x) {
  printf("\nIniciado FD");
  double ans, h = 0.00390625;
  ans = ((F(x + h) - F(x - h)) / (2 * h));
  return (ans);
}
//Segunda derivada central
double FD2(double x) {
    printf("\nIniciado FD2");
    double ans = 0.0, h = 0.00390625;
    ans = ((F(x + h)) - 2 * F(x) + F(x - h)) / MyPow(h, 2);
    return(ans);
}
//funcion que mapea la exprecion
void Tabla(float min, float max, int n) {
  float h;
  int i;
  h = (max - min) / n;
  printf("\n\t   x   \t| f(x)");
  printf("\n\t---------------------------------------");
  for (i = 0; i < n; i++) {
    printf("\n\t %f \t| %f", min + h * i, F(min + h * i));
  }
}
//funcion que mapea biseccion
void TablaTarea(double xi, double xu, int n) {
  int i;
  double xr = 0.0;
  double ea = 50;
  double es;
  double xRant;
  printf("\n\t # \t|   xi   \t|   xu   \t|   xr  \t|  F(xi)   \t|  F(xu)   \t|  F(xr)   \t|  Ea");
  printf("\n\t-----------------------------------------------------------------------------------------------------------");

  for (i = 0; i < n; i++) {
    es = Scarb(i);
    do {
      xRant = xr;
      xr = (xi + xu) / 2.0;
      if (F(xi) * F(xr) < 0.0) {
        xu = xr;
      }
      if (F(xu) * F(xr) < 0) {
        xi = xr;
      }
      ea = ErrorA(xr, xRant);
    } while (ea > es);
    printf("\n\t #%d \t|  %f  \t| %f \t| %f  \t| %f   \t|  %f   \t|  %f   \t|  %f", i, xi, xu, xr, F((float) xi), F((float) xu), F((float) xr), ea);
  }
}
//funcion de biseccion
double Biseccion(double xi, double xu, int n) {
  double xr = 0.0, xrant, es, ea = 50.0;
  es = Scarb(n);
  int i = 0;
  if (F(xi) * F(xu) > 0.0)
    printf("\n\tNo existe raiz en el intervalo");
  else
    do {
      xrant = xr;
      xr = (xi + xu) / 2.0;
      if (F(xi) * F(xr) < 0.0)
        xu = xr;
      if (F(xu) * F(xr) < 0.0)
        xi = xr;
      ea = ErrorA(xr, xrant);
      i++;
    } while (ea > es);
    printf("numero de iteraciones: %d", i);
  return (xr);
}
//metodo de newto-rapson
double NewtonRapson(double x0,int n){
    double ea = 50.0, es, x1;
    es = Scarb(n);
    printf("\n\t\t F(%f): %f",x0, F(x0));
    printf("\n\t\t FD(%f): %f",x0, FD(x0));
    printf("\n\t\tScarb: %f" , es);
    do{
        x1 = x0 - F(x0) / FD(x0);
        printf("\n\t\t En Cilco x1: %f", x1);
        printf("\n\t\t F(%f): %f",x0, F(x0));
        printf("\n\t\t FD(%f): %f",x0, FD(x0));
        ea = ErrorA(x1,x0);
        //printf("\n\t\tError Actual: %f", ea);
        x0 = x1;
    }while(ea > es);
    printf("\n\t\tError Actual: %f", ea);
    printf("\n\t\t FD(x0): %f", FD(x0));
    return(x1);
}

double NewtonRapsonMod(double x0, int n) {
    printf("Escoje el valor de P: \n\t");
    scanf("%lf",&GlobalVariable);
    printf("\nNewton Rapson Modificado inicializado");
    double ea = 50.0, es, x1;
    es = Scarb(n);
    //printf("\n\t\tF(%f)= %f", x0, F(x0));
    //printf("\n\t\tFD(%f)= %f", x0, FD(x0));
    //printf("\n\t\tFD2(%f)= %f", x0, FD2(x0));
    do{
        x1 = x0 - ((F(x0) * FD(x0)) / (MyPow(FD(x0), 2) - F(x0) * FD2(x0)));
      //  printf("\n\tx1=%f", x1);
        ea = ErrorA(x1,x0);
        x0 = x1;
        printf("\n\tEa= %f", ea);
    }while(ea > es);
    return(x1);
}

double FalsaPosicion(float xi, float xu, int n) {
  double ea = 50, es, xr = 0.0, xrant;
  int i = 0;
  es = Scarb(n);
  if (F(xi) * F(xu) > 0.0)
    printf("No existe raiz en el intervalo");
  else
    do {
    xrant = xr;
    xr = xu - (F(xu) * (xi - xu)) / (F(xi) - F(xu)); //cambia formula para xr
    if (F(xi) * F(xr) < 0.0)
      xu = xr;
    if (F(xu) * F(xr) < 0.0)
      xi = xr;
    ea = ErrorA(xr, xrant);
    i++;
  } while (ea > es);
  printf("numero de iteraciones: %d", i);
  return (xr);
}
