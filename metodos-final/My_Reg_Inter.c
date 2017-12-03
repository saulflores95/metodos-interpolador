#include "My_Reg_Inter.h"

float RegTest(float number) {
    printf("\t Number = %f ", number);
}

matriz ScanTable(void) {
  matriz a;
  int r;
  a.col = 2;
  printf("\n\t Caputra una tabla de la forma x, y ej.  2.0, 2.3");
  printf("\n\t Caputura una matriz de maximo %d pares ordenados \n\n", RMAX);
  do {
    printf("\t Numero de renglones: ");
    scanf("%d", & a.ren);
  } while (a.ren > RMAX || a.ren <= 0);

  printf("\n\n");
  for (r = 0; r < a.ren; r++) {
    printf("\tx%d, y%d= ", r + 1, r + 1);
    scanf("%f,%f", &a.mtx[r][0], &a.mtx[r][1]);
  }
  return (a);
}

matriz SumTablex(matriz t) {
  matriz sx = {4,4};
  int n;
  // n   sx  sx2 sx3
  // 00  01  02  03
  // sx  sx2 sx3 sx4
  // 10  11  12  13
  // sx2 sx3 sx4 sx5
  // 20  21  22  23
  // sx3 sx4 sx5 sx6
  // 30  31  32  33
  sx.mtx[0][0] = t.ren; //n
  for(n = 0; n < t.ren; n++){
    sx.mtx[0][1] += t.mtx[n][0]; //sumatoria en x (sx)
    sx.mtx[0][2] += MyPow(t.mtx[n][0], 2);//sx2
    sx.mtx[0][3] += MyPow(t.mtx[n][0], 3);//sx3
    sx.mtx[1][3] += MyPow(t.mtx[n][0], 4);//sx4
    sx.mtx[2][3] += MyPow(t.mtx[n][0], 5);//sx5
    sx.mtx[3][3] += MyPow(t.mtx[n][0], 6);//sx6
  }
  sx.mtx[1][0] = sx.mtx[0][1];
  sx.mtx[2][0] = sx.mtx[1][1] = sx.mtx[0][2];
  sx.mtx[3][0] = sx.mtx[1][2] = sx.mtx[2][1] = sx.mtx[0][3];
  sx.mtx[3][1] = sx.mtx[2][1] = sx.mtx[1][3];
  sx.mtx[3][2] = sx.mtx[2][3];
  return(sx);
}
  //suma de xy
matriz SumTablexy(matriz t){
  matriz sxy = { 4, 1 };
  int n;
  //sy   00
  //sxy  10
  //sx2y  20
  //sx3y 30
  for(n = 0; n < t.ren; n++) {
    sxy.mtx[0][0] += t.mtx[n][1]; //sy
    sxy.mtx[1][0] += t.mtx[n][0] * t.mtx[n][1]; //sxy
    sxy.mtx[2][0] += MyPow(t.mtx[n][0], 2) * t.mtx[n][1];//sx2y
    sxy.mtx[3][0] += MyPow(t.mtx[n][0], 3) * t.mtx[n][1];//sx3y
  }
  return(sxy);
}
  //Regresion
matriz Regresion(matriz t, int n) {
  matriz a = {0, 0}, sx, sxy;
  if (n <= 0 || n > 3) {
    printf("\n\tNo se puede realizar la regresion min 1 max 3");
  } else {
    n++;
    sx = SumTablex(t);
    sxy = SumTablexy(t);
    sx.ren = sxy.ren = sx.col = n;
    a = Cramer(sx, sxy);
  }
  return (a);
}
  //Regresion de Lineal
float RegLin(matriz t, float x) {
    float ans = 0.0;
    matriz a;
    a = Regresion(t, 1);
    printf("\n\tLa regresion lineal es:");
    if (a.mtx[1][0] >= 0.0)
        printf("\n\t y = %f + %f x", a.mtx[0][0], a.mtx[1][0]);
    else
        printf("\n\t y = %f - %f x", a.mtx[0][0], -a.mtx[1][0]);
    ans = a.mtx[0][0] + a.mtx[1][0] * x;
    printf("\n\t La evaluacion de Y(%f) = %f", x, ans);
    return (ans);
}
  //Regresion de Potencia
float RegPow(matriz t, float x){
  float ans = 0.0;
  matriz a;
  int signo = 0, n;
  for(n = 0; n < t.ren; n++) {
    if(t.mtx[n][0] < 0.0 || t.mtx[n][1] < 0.0) {
      signo = 1;
      break;
    }
    else {
      t.mtx[n][0] = log(t.mtx[n][0]);
      t.mtx[n][1] = log(t.mtx[n][1]);
    }
  }
  if(signo == 0) {
    a = Regresion(t, 1);
    printf("\n\t La regresion potencia es: ");
    printf("\n\t y = %f * X ^ %f", log(a.mtx[0][0]), a.mtx[1][0]);
    ans = a.mtx[0][0] * pow(x, a.mtx[1][0]);
    printf("\n\t La evaluacion de Y(%f) = %f", x, ans);
  }
  return(ans);
}
    //Regresion Lagrange
float Lagrang(matriz txy, float x) {
  int i, j;
  int ren = txy.ren;
  float ans = 0.0, L;

  for(i = 0; i < ren - 1; i++) {
    L = 1.0;
    for(j = 0; j < ren - 1; j++) {
      if(j != i)
        L = L * (x - txy.mtx[j][0])/(txy.mtx[i][0] - txy.mtx[j][0]);
    }
    ans = ans + L * (txy.mtx[i][1]);

  }
  return(ans);
}
    //Regresion Cuadratica
float RegCuad(matriz t, float x) {
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
    return (ans);
  }
    //Regresion Cubica
float RegCubica(matriz t, float x) {
    float ans = 0.0;
    matriz a;
    a = Regresion(t, 3);
    PrintMtx(a);
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
    return (ans);
  }
    //Regresion Exponencial
float RegExp(matriz t, float x) {
  matriz a;
  float ans = 0.0;
  int n;

  for (n = 0; n < t.ren; n++)
    t.mtx[n][1] = LogNat(t.mtx[n][1], 30);

  a = Regresion(t, 1);
  a.mtx[0][0] = Expo(a.mtx[0][0], 30);
  printf("\n\tLa regresion exponencial es: ");
  printf("\n\t y = %f * e^%fx", a.mtx[0][0], a.mtx[1][0]);
  ans = a.mtx[0][0] * Expo(a.mtx[1][0], 30);
  printf("\n\t La evaluacion de y(%.3f) = %f\n", x, ans);
  return (ans);
}
    //Regresion Lograritimica
float RegLog(matriz t, float x) {
  matriz a;
  float ans = 0.0;
  int n;

  for (n = 0; n < t.ren; n++)
    t.mtx[n][0] = LogNat(t.mtx[n][0], 20);

  a = Regresion(t, 1);
  printf("\n\tLa regresion logaritmica es: ");
  printf("\n\t y = %f * %f*Ln(x)", a.mtx[0][0], a.mtx[1][0]);
  ans = a.mtx[0][0] * a.mtx[1][0] * LogNat(x, 20);
  printf("\n\t La evaluacion de y(%.3f) = %f\n", x, ans);

  return (ans);

}
    //Regresion Inversa
float RegInv(matriz t, float x) {
  matriz a;
  float ans = 0.0;
  int n;

  for (n = 0; n < t.ren; n++)
    t.mtx[n][0] = (1 / t.mtx[n][0]);

  a = Regresion(t, 1);
  printf("\n\tLa regresion inversa es: ");
  printf("\n\t y = %f + %f/x", a.mtx[0][0], a.mtx[1][0]);
  ans = a.mtx[0][0] + a.mtx[1][0] / x;
  printf("\n\t La evaluacion de y(%.3f) = %f\n", x, ans);
  return (ans);
}

void PrintTable(matriz a) {
  PrintMtx(a);
}
