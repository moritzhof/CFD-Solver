#include"helper.h"
#include"boundary_val.h"

void boundaryvalues(int imax, int jmax, double **U, double **V){
  for(int j=1; j<=jmax; ++j){
    //Equation 15
    U[0][j] = 0.0;
    U[imax][j] = 0.0;
    // Equation 16
    V[0][j] = -V[1][j];
    V[imax+1][j] = -V[imax][j];

  }

  for(int i = 0; i<=imax; ++i){
    //equation 15
    V[i][0]=0.0;
    V[i][jmax] = 0.0;

    //equation 16
    U[i][0] = -U[i][1];
    U[i][jmax+1] = 2.0 - U[i][jmax]; // moving wall---here would be 2 times U_wall.
  }
}//end of function.
