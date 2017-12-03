#ifndef MY_MATRIX_H_INCLUDED
#define MY_MATRIX_H_INCLUDED

#include "My_Math.h"
#define RMAX 20
#define CMAX 10

typedef struct {
  int ren, col;
  float mtx[RMAX][CMAX];
} matriz;

//Funcion que lee una matriz
matriz ScanfMtx(void);
//funcion que imprime una matriz
void PrintMtx(matriz);
//matriz de ceros
matriz ZerosMtx(int, int);
//determinante de una matriz
double Determinante(matriz);
//Suma de Matriz
matriz SumaMtz(matriz, matriz);
//Resta Matriz
matriz RestaMtz(matriz, matriz);
//Multiplicacion de Matriz
matriz MultMtx(matriz, matriz);
//Verifica MTX
int VerificaMtx(matriz a, matriz b);
//Regla de cramer
matriz Cramer(matriz, matriz);
//descomposicion de Gauss
matriz DesGauss(matriz, matriz);
//Gauss - Jordan
matriz GaussJordan(matriz, matriz);
//Matriz indentidas
matriz EyeMtx(int);
//Inversa de una matris por Gauss-Jordan
matriz InvMtx(matriz);
//Metodo iterativo para sistemas de ecuaciones
matriz IteraMtx(matriz, matriz, matriz, int, int);
//Solucion de un sistema de ecuaciones por Jacobi
matriz JacobiMtx(matriz, matriz, matriz, int);
//Solucion de un sistema de ecuaciones por Gauss-Seidel
matriz GaussSeidelMtx(matriz, matriz, matriz, int);
//Transpuesta de una matriz
matriz TranspuestaMtx(matriz);
//Cambio de renglon y columna
matriz SwitchMtx(matriz a, int flag);
#endif // MY_MATRIX_H_INCLUDED
