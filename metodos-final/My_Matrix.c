#include "My_Matrix.h"
  //Funcion que lee una matriz
matriz ScanfMtx(void) {
  matriz a;
  int r, c;
  printf("\n\t Caputura una matriz de maximo %d x %d \n\n", RMAX, CMAX);
  do {
    printf("\t Numero de renglones: ");
    scanf("%d", & a.ren);
  } while (a.ren > RMAX || a.ren <= 0);

  do {
    printf("\t Numero de columnas: ");
    scanf("%d", & a.col);
  } while (a.col > CMAX || a.col <= 0);

  printf("\n\n");
  for (r = 0; r < a.ren; r++)
    for (c = 0; c < a.col; c++) {
      printf("\ta[%d][%d]", r + 1, c + 1);
      scanf("%f", & a.mtx[r][c]);
    }
  return (a);
}
  //funcion que imprime una matriz
void PrintMtx(matriz a) {
  int r, c;
  if (a.col > CMAX || a.col <= 0)
    a.col = 0;
  if (a.ren > RMAX || a.ren <= 0)
    a.ren = 0;
  printf("\n");
  for (r = 0; r < a.ren; r++) {
    for (c = 0; c < a.col; c++)
      printf("\t%.5f ", a.mtx[r][c]);
    printf("\n");
  }
}
  //matriz de ceros
matriz ZerosMtx(int r, int c) {
  matriz a = {
    r,
    c
  };
  if (a.col > CMAX || a.col <= 0)
    a.col = 0;
  if (a.ren > RMAX || a.ren <= 0)
    a.ren = 0;
  return (a);
}
  //Determinante de matriz
double Determinante(matriz a) {
  matriz maux;
  double d = 0.0;
  int r, k, c, n;
  if (a.ren != a.col)
    printf("\n\tError... La matriz debe ser cuadrada");
  else {
    if (a.ren == 1)
      d = a.mtx[0][0];
    else {
      for (n = 0; n < a.ren; n++) {
        maux = ZerosMtx(a.ren - 1, a.col - 1);
        r = 0;
        for (k = 0; k < a.ren; k++)
          if (n != k) {
            for (c = 1; c < a.col; c++)
              maux.mtx[r][c - 1] = a.mtx[k][c];
            r++;
          }
        //Esto se va a quitar
          //PrintMtx(maux);
        //Quitar antes de esto
        d += a.mtx[n][0] * MyPow(-1.0, n) * Determinante(maux);
      }
    } //a. ren == 1
  } //If a.ren != a.col
  return (d);
}
  //Suma de Matriz
matriz SumaMtz(matriz a, matriz b) {
  int i, j;
  if (a.col > CMAX || a.col <= 0)
    a.col = 0;
  if (a.ren > RMAX || a.ren <= 0)
    a.ren = 0;
  for (i = 0; i < a.ren; i++) {
    for (j = 0; j < a.col; j++)
      a.mtx[i][j] = a.mtx[i][j] + b.mtx[i][j];
  }
  return a;
}
  //Resta de Matriz
matriz RestaMtz(matriz a, matriz b) {
  int i, j;
  if (a.col > CMAX || a.col <= 0)
    a.col = 0;
  if (a.ren > RMAX || a.ren <= 0)
    a.ren = 0;
  for (i = 0; i < a.ren; i++) {
    for (j = 0; j < a.col; j++)
      a.mtx[i][j] = a.mtx[i][j] - b.mtx[i][j];
  }
  return(a);
}
  //Multiplicacion de matriz
matriz MultMtx(matriz a, matriz b) {
  if (a.col > CMAX || a.col <= 0)
    a.col = 0;
  if (a.ren > RMAX || a.ren <= 0)
    a.ren = 0;
  if(a.col == b.ren) {
    int i, j;
    for (i = 0; i < a.ren; i++) {
      for (j = 0; j < a.col; j++)
        a.mtx[i][j] = a.mtx[i][j] * b.mtx[i][j];
    }
    return(a);
  }else {
    printf("\n\t Matriz no compatible");
  }
}
  //Verificacion de matriz
int VerificaMtx(matriz a, matriz b) {
  int ans = 0;
  if (a.ren != a.col) {
    printf("\n\t La Matriz A no es cuadrada");
    ans = 1;
  }
  if (a.ren != b.ren) {
    printf("\n\t La matriz A y B tiene n numero distinto de renglones");
    ans = 2;
  }
  if (1 - (b.col == 1 || b.col == a.col)) {
    printf("\n\t Las columnas de B deben ser 1 o de la misma dimension de A");
    ans = 3;
  }
  if (Determinante(a) == 0) {
    printf("\n\t Matriz Singular");
    ans = 4;
  }
  return (ans);
}
  //Regla de Cramer
matriz Cramer(matriz a, matriz b) {
  matriz x = { 1,1 }, maux;
  int r, c;
  float da;
  int count = 1;
  if (VerificaMtx(a, b) != 0)
    printf("\n\t Nose puede resolver regla de cramer");
  else {
    da = Determinante(a);
    x = ZerosMtx(b.ren, 1);
    for (c = 0; c < a.col; c++) {
      maux = a;
      for (r = 0; r < a.ren; r++)
        maux.mtx[r][c] = b.mtx[r][0];
      x.mtx[c][0] = Determinante(maux) / da;
    }
  }
  return (x);
}
  //descomposicion de Gauss
matriz DesGauss(matriz a, matriz b) {
  int k, r, c, i;
  if (VerificaMtx(a, b) != 0) {
    printf("\n\t Nose puede resolver por Gauss");
    b = ZerosMtx(1, 1);
  } else {
    for (k = 0; k < a.ren - 1; k++){
      for (r = k + 1; r < a.ren; r++) {
        b.mtx[r][0] -= a.mtx[r][k] * b.mtx[k][0] / a.mtx[k][k];
        for (c = a.col - 1; c >= k; c--)
          a.mtx[r][c] -= a.mtx[r][k] * a.mtx[k][c] / a.mtx[k][k];
      }
    }
    printf("Ans \n");
    PrintMtx(a);
    printf("------------------\n");

    for (i = a.ren - 1; i >= 0; i--) {
      for (k = i + 1; k < a.ren; k++)
        b.mtx[i][0] -= a.mtx[i][k] * b.mtx[k][0];
      b.mtx[i][0] /= a.mtx[i][i];
    }
  }
  return (b);
}
  //Gauss - Jordan
matriz GaussJordan(matriz a, matriz b) {
  int k, r, c;
  if (VerificaMtx(a, b) != 0) {
    printf("\n\t Nose puede resolver por Gauss - Jordan");
    b = ZerosMtx(1, 1);
  } else {
    for (k = 0; k < a.ren; k++) {
      for (c = 0; c < b.col; c++)
        b.mtx[k][c] /= a.mtx[k][k];
      for (c = a.col - 1; c >= k; c--)
        a.mtx[k][c] /= a.mtx[k][k];
      for (r = 0; r < a.ren; r++)
        if (r != k) {
          for (c = 0; c < b.col; c++)
            b.mtx[r][c] -= a.mtx[r][k] * b.mtx[k][c];
          for (c = a.col - 1; c >= k; c--)
            b.mtx[r][c] -= a.mtx[r][c] * a.mtx[k][c];
        } //Fin del if
    } //Fin del for de variable k
  } //Fin de if else
  return (b);
}
  //Matriz indentidad
matriz EyeMtx(int n) {
  matriz a = {
    n,
    n
  };
  if (a.ren > CMAX || a.ren <= 0){
    a.col = 0;
    a.ren = 0;
  }
  for(n = 0; n < a.ren; n++)
    a.mtx[n][n] = 1.0;
  return (a);
}
  //Inversa de una matris por Gauss-Jordan
matriz InvMtx(matriz a) {
    matriz ide, ans;
    ide = EyeMtx(a.ren);
    if (VerificaMtx(a, ide) != 0) {
      printf("\n\tNo se puede obtener la inversa de la matriz");
      ans = ZerosMtx(1, 1);
    } else
      ans = GaussJordan(a, ide);
    return (ans);
  }
  //Metodo iterativo para sistemas de ecuaciones
matriz IteraMtx(matriz a, matriz b, matriz x0, int n, int flag) {
    matriz x1;
    float es, ea = 50.0;
    int i, j, k;
    int count = 0;
    es = Scarb(n);
    if (VerificaMtx(a, b) != 0 || VerificaMtx(a, x0) != 0) {
      printf("\n\tNo se puede resolver la matriz");
      x1 = ZerosMtx(1, 1);
    } else {
      while (ea > es) {
        x1 = b;
        for (i = 0; i < a.ren; i++) {
          for (k = 0; k <= i - 1; k++)
            if (flag == 0)
              x1.mtx[i][0] -= a.mtx[i][k] * x0.mtx[k][0]; //Jacobi
            else
              x1.mtx[i][0] -= a.mtx[i][k] * x1.mtx[k][0]; //Gauss-Seidel
          for (j = i + 1; j < a.ren; j++)
            x1.mtx[i][0] -= a.mtx[i][j] * x0.mtx[j][0];
          x1.mtx[i][0] /= a.mtx[i][i];
          count++;
        } //fin de i
        ea = ErrorA(x1.mtx[0][0], x0.mtx[0][0]);
        x0 = x1;
      } //fin de while
      printf("\n\t Cantidad de iteraciones %d \n", count);
    }
    return (x1);
  }
  //Solucion de un sistema de ecuaciones por Jacobi
matriz JacobiMtx(matriz a, matriz b, matriz x0, int n) {
    return (IteraMtx(a, b, x0, n, 0));
  }
  //Solucion de un sistema de ecuaciones por Gauss-Seidel
matriz GaussSeidelMtx(matriz a, matriz b, matriz x0, int n) {
  return (IteraMtx(a, b, x0, n, 1));
}
  //Transpuesta de una matriz
matriz TranspuestaMtx(matriz a) {
  int i, j;
  matriz ans;
  ans.ren = a.col;
  ans.col = a.ren;
  for(i = 0; i < a.ren; i++)
    for(j = 0; j < a.ren; j++){
      ans.mtx[j][i] = a.mtx[i][j];
    }
  return(ans);
}

matriz SwitchMtx(matriz a, int flag) {
  int i;
  int swp1, swp2;
  matriz newMatrix = a;
  if(flag == 1) {
    printf("\Cambia Renglon: "); scanf("%d", &swp1) ;  // Index de primer renglon
    printf("Con el Renglon: "); scanf("%d", &swp2);     // Index de segundo renglon
    for(i = 0; i < a.col; i++) {
      newMatrix.mtx[swp2 - 1][i] = a.mtx[swp1 - 1][i];
      newMatrix.mtx[swp1 - 1][i] = a.mtx[swp2 - 1][i];
    }
  }else {
    printf("\Cambia columna: "); scanf("%d", &swp1) ;  // Index de primer columna
    printf("con columna: "); scanf("%d", &swp2);     // Index de segundo columna
    for(i = 0; i < a.col; i++) {
      newMatrix.mtx[i][swp2 - 1] = a.mtx[i][swp1 - 1];
      newMatrix.mtx[i][swp1 - 1] = a.mtx[i][swp2 - 1];
    }
  }
  return(newMatrix);
}

void crout(matriz a, matriz b, matriz c, int n) {
  int i, j, k;
  double sum = 0;

  for (i = 0; i < n; i++) {
    c.mtx[i][i] = 1;
  }

  for (j = 0; j < n; j++) {
    for (i = j; i < n; i++) {
      sum = 0;
      for (k = 0; k < j; k++) {
        sum = sum + b.mtx[i][k] * c.mtx[k][j];
      }
      b.mtx[i][j] = a.mtx[i][j] - sum;
    }

    for (i = j; i < n; i++) {
      sum = 0;
      for (k = 0; k < j; k++) {
        sum = sum + b.mtx[j][k] * c.mtx[k][i];
      }
      if (b.mtx[j][j] == 0) {
        printf("det(b) close to 0!\n Can't divide by 0...\n");
        exit(EXIT_FAILURE);
      }
      c.mtx[j][i] = (a.mtx[j][i] - sum) / b.mtx[j][j];
    }
  }
  PrintMtx(c);
}

