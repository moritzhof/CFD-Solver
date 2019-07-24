#include "boundary_val.h"
#include "helper.h"
#include <stdio.h>



void boundaryvalues(int imax, int jmax, double **U, double **V){
  for(int j=1; j<=jmax; ++j){
    //Equation 15
    U[0][j] = 0.0;
    U[imax][j] = 0.0;
    // Equation 16
    V[0][j] = -V[1][j];
    V[imax+1][j] = -V[imax][j];

  }

  for(int i = 1; i<=imax; ++i){
    //equation 15
    V[i][0]= 0;
    V[i][jmax] = 0.0;

    //equation 16
    U[i][0] = -U[i][1];
    U[i][jmax+1] = 2 - U[i][jmax]; // moving wall---here would be 2 times U_wall.
  }
}

void boundaryvalues_sw(
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  int **flag
)
{
  //Obstacle boundary Conditions
    // for(int i = 1; i <= imax; ++i){
    //      MASS[i][0] = MASS[i][1];
    //     XMOM[i][0] =  XMOM[i][1];
    //     YMOM[i][0] = -YMOM[i][1];
    //     // MASS[i][70] =  MASS[i][71];
    //     // XMOM[i][70] =  XMOM[i][71];
    //     // YMOM[i][70] = -YMOM[i][71];
    //
    //     // MASS[i][71] =  MASS[i][70];
    //     // XMOM[i][71] =  XMOM[i][70];
    //     // YMOM[i][71] = -YMOM[i][70];
    //
    //     MASS[i][jmax+1] =  MASS[i][jmax];
    //     XMOM[i][jmax+1] =  XMOM[i][jmax];
    //     YMOM[i][jmax+1] = -YMOM[i][jmax];
    // }
    //
    // //Boundary Condition
    //     for (int j = 1;j <= jmax; ++j){
    //       MASS[0][j] = MASS[1][j];
    //       YMOM[0][j] = YMOM[1][j];
    //       XMOM[0][j] = -XMOM[1][j];
    //
    //       MASS[imax+1][j] = MASS[imax][j];
    //       YMOM[imax+1][j] = YMOM[imax][j];
    //        XMOM[imax+1][j] = -XMOM[imax][j];
    //     }


        for(int i = 0; i <= imax+1; ++i){
          for (int j = 0;j <= jmax+1; ++j){

            if(flag[i][j] & B_O) {
              XMOM[i][j] = -XMOM[i+1][j];
              YMOM[i][j] = YMOM[i+1][j];
              MASS[i][j] = MASS[i+1][j];

            }
            if(flag[i][j] & B_W) {
              MASS[i][j] = MASS[i-1][j];
              YMOM[i][j] = YMOM[i-1][j];
              XMOM[i][j] = -XMOM[i-1][j];

            }
            if(flag[i][j] & B_N) {
              YMOM[i][j] = -YMOM[i][j+1];
              XMOM[i][j] = XMOM[i][j+1];
              MASS[i][j] = MASS[i][j+1];

            }
            if(flag[i][j] & B_S) {
              MASS[i][j] = MASS[i][j-1];
              XMOM[i][j] = XMOM[i][j-1];
              YMOM[i][j] = -YMOM[i][j-1];

            }
          }
        }
}


void boundaryvalues_nozzle(
  int imax,
  int jmax,
  double rho,
  double ux,
  double uy,
  double p,
  double k,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  int **flag
)
{


  double U,U1,U2,RHO,RHO1,RHO2,V,V1,V2;
  double C,M;

  double p_exit = p;


//Outflow
  for(int j = 1;j <= jmax; ++j){

    //Extrapolate Rho from inside
    RHO1 = MASS[imax-1][j];
    RHO2 = MASS[imax][j];

    RHO = 2*RHO2 - RHO1;

    //Extrapolate U & V
    U1 = XMOM[imax-1][j]/MASS[imax-1][j];
    U2 = XMOM[imax][j]/MASS[imax][j];

    V1 = YMOM[imax-1][j]/MASS[imax-1][j];
    V2 = YMOM[imax][j]/MASS[imax][j];

    U = 2*U2 - U1;
    V = 2*V2 - V1;

    C =sqrt(k*p_exit/RHO);
    M = sqrt(U*U + V*V)/C;



 if(M <=1){

    MASS[imax+1][j] = RHO;
    XMOM[imax+1][j] = U*RHO;
    YMOM[imax+1][j] = V*RHO;
    ENER[imax+1][j] = 0.5*RHO*(U*U+V*V)+p_exit/(k-1);
  }
  }

  for(int i = 0; i <= imax+1; ++i){
    for (int j = 0;j <= jmax+1; ++j){
      // if(flag[i][j] & B_O) XMOM[i][j] = -XMOM[i+1][j];
      // if(flag[i][j] & B_W) XMOM[i][j] = -XMOM[i-1][j];
      if(flag[i][j] & B_O) {
        XMOM[i][j] = -XMOM[i+1][j];
        // YMOM[i][j] = 0;
      }
      if(flag[i][j] & B_W) {
        XMOM[i][j] = -XMOM[i-1][j];
        // YMOM[i][j] = 0;
      }
      if(flag[i][j] & B_N) {
        YMOM[i][j] = -YMOM[i][j+1];
        // XMOM[i][j] = 0;
      }
      if(flag[i][j] & B_S) {
        YMOM[i][j] = -YMOM[i][j-1];
        // XMOM[i][j] = 0;
      }
    }
  }
}
