#include "My_Math.h"
#include "My_Roots.h"
#include "My_Poli_Roots.h"
#include "My_Matrix.h"
#include "My_Reg_Inter.h"
#include "My_Final_Project.h"
int main() {
    //sim1/3 -> fi -> normal
    //f -> sim1/3 -p
    //nr f
    //d/d=f
    //Proyecto Final Metodos Numericos
    int option;
    double value;
    int exit = 1;
    do {
        printf("***********************************************************\n");
        printf("**********Bienvenido al Proyecto Final de metodos**********\n");
        printf("***********************************************************\n");
        printf("%f", ElPow(2, 3));
        printf("\n\t Seleccione opcion: ");
        printf("\n\t 1. Best Fit ");
        printf("\n\t 2. Integrador \n\t");
        scanf("%d", &option);
        if(option == 1) {
            printf("\n\tInterpolador con Best Fitt\n\t");
            //Interpolador(0.0, 2.0, 1.0, 4.0);
            //Interpolador(1, -2, 1.0, 10);
            Interpolador(0, 1, 1.0, 20);
        }
        if(option == 2) {
            printf("\n\tIntegrador 2000\n");
            Integrador(0.5, 0.5);
            // printf("\n\tNewton: %f", NewtonRapson(0.5, 5));
            //printf("\n\tNewton - MOD: %f", NewtonRapsonMod(0.5, 5));
        }
        printf("\n\tQuieres continuar?");
        printf("\n\t 1. Continuar");
        printf("\n\t 2. Salir \n\t");
        scanf("%d", &exit);

    }while(exit == 1);
    return 0;
}
