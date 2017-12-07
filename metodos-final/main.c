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
    printf("\tInterpolador con Best Fitt\n");
    //Interpolador(0.0, 2.0, 1.0, 4.0);
    //Interpolador(1, -2, 1.0, 10);
//    printf("\n\tNewton: %f", NewtonRapson(0.5, 5));
    printf("\n\tNewton - MOD: %f", NewtonRapsonMod(0.5, 5));
    return 0;
}
