#include "uvp.h"
#include <stdio.h>
#include "helper.h"
#include <math.h>


void x_llf_ssc(
  double* dt,
  int imax,
  int jmax,
  double rho,
  double ux,
  double uy,
  double p,
  double k,
  double dx,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  double **resMass,
  double **resXMom,
  double **resYMom,
  double **resEner,
  int **flag
)
{
  //Free stream boundary conditions
  for(int j = 1;j <= jmax; ++j){
    MASS[0][j] = rho;
    XMOM[0][j] = ux*rho;
    YMOM[0][j] = uy*rho;
    ENER[0][j] = 0.5*rho*(ux*ux+uy*uy)+p/(k-1);

    MASS[imax+1][j] = rho;
    XMOM[imax+1][j] = ux*rho;
    YMOM[imax+1][j] = uy*rho;
    ENER[imax+1][j] = 0.5*rho*(ux*ux+uy*uy)+p/(k-1);
  }
//Obstacle boundary Conditions
  for(int i = 1; i <= imax; ++i){
    for (int j = 1;j <= jmax; ++j){
      if(flag[i][j] & B_O) {
        XMOM[i][j] = -XMOM[i+1][j];
        MASS[i][j] = MASS[i+1][j];
        YMOM[i][j] = YMOM[i+1][j];
        ENER[i][j] = ENER[i+1][j];
      }
      if(flag[i][j] & B_W) {
        MASS[i][j] = MASS[i-1][j];
        XMOM[i][j] = -XMOM[i-1][j];
        YMOM[i][j] = YMOM[i-1][j];
        ENER[i][j] = ENER[i-1][j];
      }
    }
  }


//Calculate x-direction flux F Right and Left
double xFluxR0,xFluxL0;
double xFluxR1,xFluxL1;
double xFluxR2,xFluxL2;
double xFluxR3,xFluxL3;
double pressR,pressL;

double flux0,flux1,flux2,flux3;

double lambdaR,lambdaL,maxLambda;
double DX = 1/dx;


for(int i = 1; i <= imax+1; ++i){
  for (int j = 1;j <= jmax+1; ++j){
    if(flag[i][j] & FLUID || flag[i][j] & B_W){
      pressL = (ENER[i-1][j] - 0.5*(XMOM[i-1][j]*XMOM[i-1][j]+YMOM[i-1][j]*YMOM[i-1][j])/MASS[i-1][j])*(k-1);


      xFluxL0 = XMOM[i-1][j];
      xFluxL1 = XMOM[i-1][j]*XMOM[i-1][j]/MASS[i-1][j] + pressL;
      xFluxL2 = XMOM[i-1][j]*YMOM[i-1][j]/MASS[i-1][j];
      xFluxL3 = XMOM[i-1][j]/MASS[i-1][j]*(k*pressL/(k-1) + 0.5*(XMOM[i-1][j]*XMOM[i-1][j] + YMOM[i-1][j]*YMOM[i-1][j])/MASS[i-1][j]);

      pressR = (ENER[i][j] - 0.5*(XMOM[i][j]*XMOM[i][j]+YMOM[i][j]*YMOM[i][j])/MASS[i][j])*(k-1);

      xFluxR0 = XMOM[i][j];
      xFluxR1 = XMOM[i][j]*XMOM[i][j]/MASS[i][j] + pressR;
      xFluxR2 = XMOM[i][j]*YMOM[i][j]/MASS[i][j];
      xFluxR3 = XMOM[i][j]/MASS[i][j]*(k*pressR/(k-1) + 0.5*(XMOM[i][j]*XMOM[i][j] + YMOM[i][j]*YMOM[i][j])/MASS[i][j]);

      lambdaL = fabs(XMOM[i-1][j]/MASS[i-1][j]) + sqrt(k*pressL/MASS[i-1][j]);
      lambdaR = fabs(XMOM[i  ][j]/MASS[i  ][j]) + sqrt(k*pressR/MASS[i  ][j]);

      maxLambda = fmax(lambdaL,lambdaR);

      flux0 =  0.5*(xFluxL0 + xFluxR0) - 0.5*maxLambda*(MASS[i][j] - MASS[i-1][j]);
      flux1 =  0.5*(xFluxL1 + xFluxR1) - 0.5*maxLambda*(XMOM[i][j] - XMOM[i-1][j]);
      flux2 =  0.5*(xFluxL2 + xFluxR2) - 0.5*maxLambda*(YMOM[i][j] - YMOM[i-1][j]);
      flux3 =  0.5*(xFluxL3 + xFluxR3) - 0.5*maxLambda*(ENER[i][j] - ENER[i-1][j]);

      resMass[i  ][j] += DX*flux0;
      resMass[i-1][j] -= DX*flux0;

      resXMom[i  ][j] += DX*flux1;
      resXMom[i-1][j] -= DX*flux1;

      resYMom[i  ][j] += DX*flux2;
      resYMom[i-1][j] -= DX*flux2;

      resEner[i  ][j] += DX*flux3;
      resEner[i-1][j] -= DX*flux3;

      *dt = fmin(*dt,cfl*dx/maxLambda);
      }
    }
  }
}

void y_llf_ssc(
  double* dt,
  int imax,
  int jmax,
  double rho,
  double ux,
  double uy,
  double p,
  double k,
  double dy,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  double **resMass,
  double **resXMom,
  double **resYMom,
  double **resEner,
  int **flag
)
{
  //Free stream boundary conditions
  for(int i = 1;i <= imax; ++i){
    MASS[i][0] = rho;
    XMOM[i][0] = ux*rho;
    YMOM[i][0] = uy*rho;
    ENER[i][0] = 0.5*rho*(ux*ux+uy*uy)+p/(k-1);

    MASS[i][jmax+1] = rho;
    XMOM[i][jmax+1] = ux*rho;
    YMOM[i][jmax+1] = uy*rho;
    ENER[i][jmax+1] = 0.5*rho*(ux*ux+uy*uy)+p/(k-1);
  }
//Obstacle boundary Conditions
  for(int i = 1; i <= imax; ++i){
    for (int j = 1;j <= jmax; ++j){
      if(flag[i][j] & B_N) {
        MASS[i][j] = MASS[i][j+1];
        XMOM[i][j] = XMOM[i][j+1];
        YMOM[i][j] = -YMOM[i][j+1];
        ENER[i][j] = ENER[i][j+1];
      }
      if(flag[i][j] & B_S) {
        MASS[i][j] = MASS[i][j-1];
        XMOM[i][j] = XMOM[i][j-1];
        YMOM[i][j] = -YMOM[i][j-1];
        ENER[i][j] = ENER[i][j-1];
      }
    }
  }


//Calculate x-direction flux F Right and Left
double yFluxT0,yFluxB0;
double yFluxT1,yFluxB1;
double yFluxT2,yFluxB2;
double yFluxT3,yFluxB3;
double pressT,pressB;

double flux0,flux1,flux2,flux3;

double lambdaT,lambdaB,maxLambda;
double DY = 1/dy;


for(int i = 1; i <= imax+1; ++i){
  for (int j = 1;j <= jmax+1; ++j){
    if(flag[i][j] & FLUID || flag[i][j] & B_S){
      pressB = (ENER[i][j-1] - 0.5*(XMOM[i][j-1]*XMOM[i][j-1]+YMOM[i][j-1]*YMOM[i][j-1])/MASS[i][j-1])*(k-1);

      yFluxB0 = YMOM[i][j-1];
      yFluxB1 = YMOM[i][j-1]*XMOM[i][j-1]/MASS[i][j-1];
      yFluxB2 = YMOM[i][j-1]*YMOM[i][j-1]/MASS[i][j-1] + pressB;
      yFluxB3 = YMOM[i][j-1]/MASS[i][j-1]*(k*pressB/(k-1) + 0.5*(XMOM[i][j-1]*XMOM[i][j-1] + YMOM[i][j-1]*YMOM[i][j-1])/MASS[i][j-1]);

      pressT = (ENER[i][j] - 0.5*(XMOM[i][j]*XMOM[i][j]+YMOM[i][j]*YMOM[i][j])/MASS[i][j])*(k-1);

      yFluxT0 = YMOM[i][j];
      yFluxT1 = YMOM[i][j]*XMOM[i][j]/MASS[i][j];
      yFluxT2 = YMOM[i][j]*YMOM[i][j]/MASS[i][j] + pressT;
      yFluxT3 = YMOM[i][j]/MASS[i][j]*(k*pressT/(k-1) + 0.5*(XMOM[i][j]*XMOM[i][j] + YMOM[i][j]*YMOM[i][j])/MASS[i][j]);

      lambdaB = fabs(YMOM[i][j-1]/MASS[i][j-1]) + sqrt(k*pressB/MASS[i][j-1]);
      lambdaT = fabs(YMOM[i][j  ]/MASS[i][j  ]) + sqrt(k*pressT/MASS[i][j  ]);

      maxLambda = fmax(lambdaB,lambdaT);

      flux0 =  0.5*(yFluxB0 + yFluxT0) - 0.5*maxLambda*(MASS[i][j] - MASS[i][j-1]);
      flux1 =  0.5*(yFluxB1 + yFluxT1) - 0.5*maxLambda*(XMOM[i][j] - XMOM[i][j-1]);
      flux2 =  0.5*(yFluxB2 + yFluxT2) - 0.5*maxLambda*(YMOM[i][j] - YMOM[i][j-1]);
      flux3 =  0.5*(yFluxB3 + yFluxT3) - 0.5*maxLambda*(ENER[i][j] - ENER[i][j-1]);

      resMass[i][j  ] += DY*flux0;
      resMass[i][j-1] -= DY*flux0;

      resXMom[i][j  ] += DY*flux1;
      resXMom[i][j-1] -= DY*flux1;

      resYMom[i][j  ] += DY*flux2;
      resYMom[i][j-1] -= DY*flux2;

      resEner[i][j  ] += DY*flux3;
      resEner[i][j-1] -= DY*flux3;

      *dt = fmin(*dt,cfl*dy/maxLambda);
      }
    }
  }
}


void time_step_ssc(
  double dt,
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  double **resMass,
  double **resXMom,
  double **resYMom,
  double **resEner
)
{

  for(int i = 1; i<= imax; ++i){
    for(int j = 1; j<= jmax; ++j){
    MASS[i][j] = MASS[i][j] + dt*resMass[i][j];
    XMOM[i][j] = XMOM[i][j] + dt*resXMom[i][j];
    YMOM[i][j] = YMOM[i][j] + dt*resYMom[i][j];
    ENER[i][j] = ENER[i][j] + dt*resEner[i][j];
  }
}
}

void calculate_fg_cds_fv(
  double dt,
  int imax,
  int jmax,
  double dx,
  double dy,
  double **U,
  double **V,
  double **F,
  double **G,
  int **flag
)
{


double DY = 1/dy;
double DX = 1/dx;
double Re = 100;
double RE = 1/Re;

double ue,uw,vn,vs;
double Ae,Aw,An,As,Ap;


for(int i = 1; i <= imax; ++i){
  for (int j = 1;j <= jmax; ++j){
    if(flag[i][j] & FLUID){

      ue = 0.5*(U[i+1][j] + U[i][j]);
      uw = 0.5*(U[i-1][j] + U[i][j]);
      vn = 0.5*(V[i][j+1] + V[i-1][j+1]);
      vs = 0.5*(V[i-1][j] + V[i][j]);


        ue = 0.5*(U[i+1][j] + U[i][j]);
        uw = 0.5*(U[i-1][j] + U[i][j]);
        vn = 0.5*(V[i][j+1] + V[i-1][j+1]);
        vs = 0.5*(V[i-1][j] + V[i][j]);

        Ae = (+0.5*ue - RE*DX)*dy;
        Aw = (-0.5*uw - RE*DX)*dy;
        An = (+0.5*vn - RE*DY)*dx;
        As = (-0.5*vs - RE*DY)*dx;
        Ap = (0.5*ue + RE*DX)*dy + (-0.5*uw+ RE*DX)*dy + (0.5*vn + RE*DY)*dx + (-0.5*vs + RE*DY)*dx;

     F[i][j] =U[i][j] - dt*(Ae*U[i+1][j] + Aw*U[i-1][j] +An*U[i][j+1] + As*U[i][j-1] + Ap*U[i][j])/(dx*dy);
      }
    }
  }

  for(int i = 1; i <= imax; ++i){
    for (int j = 1;j <= jmax; ++j){
      if(flag[i][j] & FLUID){

        ue = 0.5*(U[i+1][j] + U[i][j]);
        uw = 0.5*(U[i-1][j] + U[i][j]);
        vn = 0.5*(V[i][j+1] + V[i-1][j+1]);
        vs = 0.5*(V[i-1][j] + V[i][j]);

        Ae = (+0.5*ue - RE*DX)*dy;
        Aw = (-0.5*uw - RE*DX)*dy;
        An = (+0.5*vn - RE*DY)*dx;
        As = (-0.5*vs - RE*DY)*dx;
        Ap = (0.5*ue + RE*DX)*dy + (-0.5*uw+ RE*DX)*dy + (0.5*vn + RE*DY)*dx + (-0.5*vs + RE*DY)*dx;

       G[i][j] = V[i][j]-dt*(Ae*V[i+1][j] +Aw*V[i-1][j] +An*V[i][j+1] + As*V[i][j-1] + Ap*V[i][j])/(dx*dy);
        }
      }
    }

}


void calculate_fg_uds_fv(
  double dt,
  int imax,
  int jmax,
  double dx,
  double dy,
  double **U,
  double **V,
  double **F,
  double **G,
  int **flag
)
{


double DY = 1/dy;
double DX = 1/dx;
double Re = 100;
double RE = 1/Re;

double ue,uw,vn,vs;
double Ae,Aw,An,As,Ap;

// double uue,uuw;
// double vus,vun;

for(int i = 1; i <= imax; ++i){
  for (int j = 1;j <= jmax; ++j){
    if(flag[i][j] & FLUID){

      ue = 0.5*(U[i+1][j] + U[i][j]);
      uw = 0.5*(U[i-1][j] + U[i][j]);
      vn = 0.5*(V[i][j+1] + V[i-1][j+1]);
      vs = 0.5*(V[i-1][j] + V[i][j]);

      Ae = (fmax( ue,0.0) - RE*DX)*dy;
      Aw = (-fmax(uw,0.0) - RE*DX)*dy;
      An = (fmax( vn,0.0) - RE*DY)*dx;
      As = (-fmax(vs,0.0) - RE*DY)*dx;
      Ap = (fmax(ue,0.0) + RE*DX)*dy + (-fmax(uw,0.0) + RE*DX)*dy + (fmax(vn,0.0) + RE*DY)*dx + (-fmax(vs,0.0)+ RE*DY)*dx;



     F[i][j] =U[i][j] -dt*(Ae*U[i+1][j] +Aw*U[i-1][j] +An*U[i][j+1] + As*U[i][j-1] + Ap*U[i][j])/(dx*dy);
      }
    }
  }

  for(int i = 1; i <= imax; ++i){
    for (int j = 1;j <= jmax; ++j){
      if(flag[i][j] & FLUID){

        ue = 0.5*(U[i+1][j] + U[i][j]);
        uw = 0.5*(U[i-1][j] + U[i][j]);
        vn = 0.5*(V[i][j+1] + V[i-1][j+1]);
        vs = 0.5*(V[i-1][j] + V[i][j]);

        Ae = (fmax( ue,0.0) - RE*DX)*dy;
        Aw = (-fmax(uw,0.0) - RE*DX)*dy;
        An = (fmax( vn,0.0) - RE*DY)*dx;
        As = (-fmax(vs,0.0) - RE*DY)*dx;
        Ap = (fmax(ue,0.0) + RE*DX)*dy + (-fmax(uw,0.0) + RE*DX)*dy + (fmax(vn,0.0) + RE*DY)*dx + (-fmax(vs,0.0)+ RE*DY)*dx;

       G[i][j] = V[i][j]-dt*(Ae*V[i+1][j] +Aw*V[i-1][j] +An*V[i][j+1] + As*V[i][j-1] + Ap*V[i][j])/(dx*dy);
        }
      }
    }
}

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
){
  double tmp, maxi1, maxi2;
  /* we rewrite dt only if tau is nonnegative, else we do nothing */
  if(tau>=0){
    tmp = 0.5*Re*(dx*dx*dy*dy)/(dx*dx+dy*dy);
    maxi1 = mmax(U, imax, jmax);
    maxi2 = mmax(V, imax, jmax);
    *dt = tau * fmin(tmp, fmin(dx/maxi1, dy/maxi2)); // (dx/mmax(U, imax, jmax), dy/mmax(V, imax, jmax));
  }
}



void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
){
  for(int i = 1; i <= imax; ++i){
    for(int j = 1; j <= jmax; ++j){
      U[i][j] = F[i][j] - (dt/dx)*(P[i][j]-P[i-1][j]);
      V[i][j] = G[i][j] - (dt/dy)*(P[i][j]-P[i][j-1]);

    }
  }

}


void x_llf_sw(
  double* dt,
  int imax,
  int jmax,
  double g,
  double dx,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **resMass,
  double **resXMom,
  double **resYMom,
  int **flag
)
{

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



//Calculate x-direction flux F Right and Left
double xFluxR0,xFluxL0;
double xFluxR1,xFluxL1;
double xFluxR2,xFluxL2;
double flux0,flux1,flux2;

double lambdaR,lambdaL,maxLambda;
double DX = 1/dx;


for(int i = 1; i <= imax+1; ++i){
  for (int j = 1;j <= jmax+1; ++j){
        // if(flag[i][j] & FLUID || flag[i][j] & B_W){
// if(flag[i][j] & OUTFLOW) g = 0;

      xFluxL0 = XMOM[i-1][j];
      xFluxL1 = XMOM[i-1][j]*XMOM[i-1][j]/MASS[i-1][j] + 0.5*g*MASS[i-1][j]*MASS[i-1][j];
      xFluxL2 = XMOM[i-1][j]*YMOM[i-1][j]/MASS[i-1][j];

      xFluxR0 = XMOM[i][j];
      xFluxR1 = XMOM[i][j]*XMOM[i][j]/MASS[i][j] +  0.5*g*MASS[i][j]*MASS[i][j];
      xFluxR2 = XMOM[i][j]*YMOM[i][j]/MASS[i][j];

      lambdaL = fabs(XMOM[i-1][j]/MASS[i-1][j]) + sqrt(g*MASS[i-1][j]);
      lambdaR = fabs(XMOM[i  ][j]/MASS[i  ][j]) + sqrt(g*MASS[i  ][j]);

      maxLambda = fmax(lambdaL,lambdaR);

      flux0 =  0.5*(xFluxL0 + xFluxR0) - 0.5*maxLambda*(MASS[i][j] - MASS[i-1][j]);
      flux1 =  0.5*(xFluxL1 + xFluxR1) - 0.5*maxLambda*(XMOM[i][j] - XMOM[i-1][j]);
      flux2 =  0.5*(xFluxL2 + xFluxR2) - 0.5*maxLambda*(YMOM[i][j] - YMOM[i-1][j]);

      resMass[i  ][j] += DX*flux0;
      resMass[i-1][j] -= DX*flux0;

      resXMom[i  ][j] += DX*flux1;
      resXMom[i-1][j] -= DX*flux1;

      resYMom[i  ][j] += DX*flux2;
      resYMom[i-1][j] -= DX*flux2;


      *dt = fmin(*dt,cfl*dx/maxLambda);
// }
    }
  }
}

void y_llf_sw(
  double* dt,
  int imax,
  int jmax,
  double g,
  double dy,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **resMass,
  double **resXMom,
  double **resYMom,
  int **flag
)
{

// //Obstacle boundary Conditions
//   for(int i = 1; i <= imax; ++i){
//        MASS[i][0] = MASS[i][1];
//       XMOM[i][0] =  XMOM[i][1];
//       YMOM[i][0] = -YMOM[i][1];
//
//       MASS[i][jmax+1] =  MASS[i][jmax];
//       XMOM[i][jmax+1] =  XMOM[i][jmax];
//       YMOM[i][jmax+1] = -YMOM[i][jmax];
//   }


//Calculate x-direction flux F Right and Left
double yFluxT0,yFluxB0;
double yFluxT1,yFluxB1;
double yFluxT2,yFluxB2;

double flux0,flux1,flux2;

double lambdaT,lambdaB,maxLambda;
double DY = 1/dy;


for(int i = 1; i <= imax+1; ++i){
  for (int j = 1;j <= jmax+1; ++j){
        // if(flag[i][j] & FLUID || flag[i][j] & B_W){
// if(flag[i][j] & OUTFLOW) g = 0;

      yFluxB0 = YMOM[i][j-1];
      yFluxB1 = YMOM[i][j-1]*XMOM[i][j-1]/MASS[i][j-1];
      yFluxB2 = YMOM[i][j-1]*YMOM[i][j-1]/MASS[i][j-1] +  0.5*g*MASS[i][j-1]*MASS[i][j-1];

      yFluxT0 = YMOM[i][j];
      yFluxT1 = YMOM[i][j]*XMOM[i][j]/MASS[i][j];
      yFluxT2 = YMOM[i][j]*YMOM[i][j]/MASS[i][j] +  0.5*g*MASS[i][j]*MASS[i][j];

      lambdaB = fabs(YMOM[i][j-1]/MASS[i][j-1]) + sqrt(g*MASS[i][j-1]);
      lambdaT = fabs(YMOM[i][j  ]/MASS[i][j  ]) + sqrt(g*MASS[i][j  ]);

      maxLambda = fmax(lambdaB,lambdaT);

      flux0 =  0.5*(yFluxB0 + yFluxT0) - 0.5*maxLambda*(MASS[i][j] - MASS[i][j-1]);
      flux1 =  0.5*(yFluxB1 + yFluxT1) - 0.5*maxLambda*(XMOM[i][j] - XMOM[i][j-1]);
      flux2 =  0.5*(yFluxB2 + yFluxT2) - 0.5*maxLambda*(YMOM[i][j] - YMOM[i][j-1]);

      resMass[i][j  ] += DY*flux0;
      resMass[i][j-1] -= DY*flux0;

      resXMom[i][j  ] += DY*flux1;
      resXMom[i][j-1] -= DY*flux1;

      resYMom[i][j  ] += DY*flux2;
      resYMom[i][j-1] -= DY*flux2;

      *dt = fmin(*dt,cfl*dy/maxLambda);
      // }
     }
   }
  }









void time_step_sw(
  double dt,
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **resMass,
  double **resXMom,
  double **resYMom,
  int **flag
)
{

  for(int i = 1; i<= imax; ++i){
    for(int j = 1; j<= jmax; ++j){
    if(flag[i][j] & FLUID){
    MASS[i][j] = MASS[i][j] + dt*resMass[i][j];
    XMOM[i][j] = XMOM[i][j] + dt*resXMom[i][j];
    YMOM[i][j] = YMOM[i][j] + dt*resYMom[i][j];
  }
  }
}
}
