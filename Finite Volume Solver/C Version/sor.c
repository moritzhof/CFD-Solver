#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  rloc = rloc/(imax*jmax);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


  /* set boundary values */
  for(i = 1; i <= imax; i++) {
    P[i][0] = P[i][1];
    P[i][jmax+1] = P[i][jmax];
  }
  for(j = 1; j <= jmax; j++) {
    P[0][j] = P[1][j];
    P[imax+1][j] = P[imax][j];
  }
}

void caluculate_p_fv(
    int imax,
    int jmax,
    int itmax,
    double dt,
    double dx,
    double dy,
    double omg,
    double eps,
    double **F,
    double **G,
    double **P,
    double **RS,
    double **resP
  )
{
    double DT = 1/dt;
    //Calculate RS
    for(int i = 1; i <= imax; ++i){
      for(int j = 1; j <= jmax; ++j){
        RS[i][j] = ((F[i+1][j] - F[i][j])*dy + (G[i][j+1] - G[i][j])*dx)*DT;

      }
    }

    double  res = 10;
    double Asp,Anp,Awp,Aep,App;
    double DX = 1/dx;
    double DY = 1/dy;
    int  it = 0;

    while (res > eps && it < itmax){
 // while (res > eps){

      res = 0;
      for(int i = 1; i <= imax; ++i){
        for(int j = 1; j <= jmax; ++j){
          //Boundary Conditions
            Asp = -dx*DY;
            if(j == 1)    Asp = 0;

            Anp = -dx*DY;
            if(j == jmax) Anp = 0;

            Awp = -dy*DX;
             if(i == 1)    Awp = 0;

            //
            Aep = -dy*DX;
            if(i == imax) Aep = 0;


            App = Awp + Aep + Asp + Anp;
            resP[i][j] = (1/App)*(App*(1-omg)*P[i][j] + omg*(Asp*resP[i][j-1]+Awp*resP[i-1][j] + Aep*P[i+1][j] + Anp*P[i][j+1] + RS[i][j]));
            res       +=  (Aep*P[i+1][j] + Awp*P[i-1][j] + Anp*P[i][j+1] + Asp*P[i][j-1] - Asp*P[i][j] + RS[i][j])*
                          (Aep*P[i+1][j] + Awp*P[i-1][j] + Anp*P[i][j+1] + Asp*P[i][j-1] - Asp*P[i][j] + RS[i][j]);
           // P[i][j] = resP[i][j];
        }
      }

      for(int i = 1; i <= imax; ++i){
        for(int j = 1; j <= jmax; ++j){
           P[i][j] = resP[i][j];
        }}

        /* set boundary values */
        for(int i = 1; i <= imax; i++) {
          P[i][0] = P[i][1];
          P[i][jmax+1] = P[i][jmax];
        }
        for(int j = 1; j <= jmax; j++) {
          P[0][j] = P[1][j];
          P[imax+1][j] = P[imax][j];
        }
      res = res/(imax*jmax);
      res = sqrt(res);
      ++it;
        // printf("res =%f",res);
      /* set boundary values */

    }
}
