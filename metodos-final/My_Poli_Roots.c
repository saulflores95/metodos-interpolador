#include "My_Poli_Roots.h"

polinomio ScanfPol(void) {
  polinomio p;
  int n;
  printf("\n\t Captura un polinomio de grado maximo %d", GMAX - 1);
  printf("\n\t De la forma a0 + a1X + a2 X^2 + ... + an X^n");
  do {
    printf("\n\t Grado del polinomio: ");
    scanf("%d", & p.grado);
  } while ((p.grado > GMAX - 1) || (p.grado < 0));
  for (n = 0; n <= p.grado; n++) {
    printf("\ta[%d] = ", n);
    scanf("%f", & p.pol);
  }
  return (p);
}

void PrintPol(polinomio p) {
  int n;
  printf("\n\t");
  if (p.pol[0] == 0)
    printf("\n\t");
  else
  if (p.pol[0] < 0)
    printf("- %f", Abs(p.pol[0]));
  else
    printf(" %f", p.pol[0]);
  for (n = 1; n <= p.grado; n++) {
    if (p.pol[n] != 0.0)
      if (n == 1) {
        if (p.pol[n] < 0)
          printf(" - %f X", Abs(p.pol[n]));
        else
          printf(" + %f X", p.pol[n]);
      } else {
        if (p.pol[n] < 0)
          printf(" - %f X^%d", Abs(p.pol[n]), n);
        else
          printf(" + %f X^%d", p.pol[n], n);
      }
  }
  printf("\n\t");
}

float PolEval(polinomio a, float x) {
  float ans = 0.0;
  int n;
  for (n = a.grado; n >= 0; n--)
    ans = ans * x + a.pol[n];
  return (ans);
}

polinomio PolDeriv(polinomio a) {
  polinomio b = {
    0
  };
  int n;
  if (a.grado > 0) {
    b.grado = a.grado - 1;
    for (n = 0; n < a.grado; n++)
      b.pol[n] = a.pol[n + 1] * (n + 1);
  }
  return (b);
}

double NewtonRPol(polinomio f, float x0, int n) {
  double x1, ea = 50.0, es;
  polinomio fd;
  es = Scarb(n);
  fd = PolDeriv(f);
  int i = 0;
  do {
    x1 = x0 - PolEval(f, x0) / PolEval(fd, x0);
    ea = ErrorA(x1, x0);
    x0 = x1;
    i++;
  } while (ea > es);
  printf("Numero de iteraciones: %d", i);
  return (x1);
}

//Division Sintetica
polinomio DivisionSint(polinomio a, float x){
  polinomio b;
  int i;
  b.grado = a.grado - 1;
  b.pol[b.grado] = a.pol[a.grado];
  for(i = b.grado - 1; i > 0; i--)
    b.pol[i] = a.pol[i + 1] + b.pol[i + 1] * x;
  return(b);
}

//Metodo de Birge-Vierta
void BirgeVieta(polinomio a) {
  int len = sizeof(a.pol) / sizeof(a.pol[0]); //Se declara len para la longitud del polinomio
  double newtonRapsonRaiz = NewtonRPol(a, -2.0, 5); //Se consigue la raiz utilizando el metodo de rapzon
  printf("\n\t la Raiz es %.5f", newtonRapsonRaiz);
  PrintPol(a);
  int i = 0;
  do {
    a = DivisionSint(a, newtonRapsonRaiz);
    printf("\tGRADO: %d", a.grado);
    PrintPol(a);
    i++;
  } while (a.grado > 0);
  printf("Numero de iteraciones: %d", i);
}
