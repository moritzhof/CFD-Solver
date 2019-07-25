#include"helper.h"
#include"init.h"

#include<iostream>
#include<string>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<iterator>
#include<algorithm>
#include<fstream>
#include<math.h>

void xLocalLaxFriedrichs(double* dt, int imax, int jmax, double rho, double ux, double uy,
                         double p, double k,  double dx, double cfl, Matrix& MASS, Matrix& XMOM,
                         Matrix& YMOM, Matrix& ENER, Matrix& resMass, Matrix& resXMom, Matrix& resYMom,
                         Matrix& resEner, intMatrix& FLAG){

    // std::cout<< "FLAG SIZE: " << FLAG.size()<< " "<< FLAG[0].size()<<std::endl;
    // std::cout<< "YMOM SIZE: " << YMOM.size()<< " "<< YMOM[0].size()<<std::endl;
    // std::cout<< "XMOM SIZE: " << XMOM.size()<< " "<< YMOM[0].size()<<std::endl;
    // std::cout<< "ENER SIZE: " << ENER.size()<< " "<< ENER[0].size()<<std::endl;
    // std::cout<< "MASS SIZE: " << MASS.size()<< " "<< MASS[0].size()<<std::endl;
    // std::cout<< "resEner SIZE: " << resEner.size()<< " "<< resEner[0].size()<<std::endl;
    // std::cout<< "resMassS: " << FLAG.size()<< " "<< FLAG[0].size()<<std::endl;
    // std::cout<< "FL

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
      if(FLAG[i][j] & State::B_O) XMOM[i][j] = -XMOM[i+1][j];
      if(FLAG[i][j] & State::B_W) XMOM[i][j] = -XMOM[i-1][j];
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
    if(FLAG[i][j] & State::FLUID || FLAG[i][j] & State::B_W){
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

      lambdaL = std::fabs(XMOM[i-1][j]/MASS[i-1][j]) + std::sqrt(k*pressL/MASS[i-1][j]);
      lambdaR = std::fabs(XMOM[i  ][j]/MASS[i  ][j]) + std::sqrt(k*pressR/MASS[i  ][j]);

      maxLambda = std::fmax(lambdaL,lambdaR);

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


void yLocalLaxFriedrichs(double* dt, int imax, int jmax, double rho, double ux, double uy,
                           double p, double k,  double dy, double cfl, Matrix& MASS, Matrix& XMOM,
                           Matrix& YMOM, Matrix& ENER, Matrix& resMass, Matrix& resXMom, Matrix& resYMom,
                           Matrix& resEner, intMatrix& FLAG){

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
                                   if(FLAG[i][j] & State::B_N) YMOM[i][j] = -YMOM[i][j+1];
                                   if(FLAG[i][j] & State::B_S) YMOM[i][j] = -YMOM[i][j-1];
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
                                 if(FLAG[i][j] & State::FLUID || FLAG[i][j] & State::B_S){
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


void timeStep(double dt, int imax, int jmax, Matrix& MASS, Matrix& XMOM, Matrix& YMOM,
              Matrix& ENER, Matrix& resMass,Matrix& resXMom, Matrix& resYMom, Matrix& resEner){
  for(int i = 1; i<=imax; ++i){
    for(int j = 1; j<=jmax; ++j){
    MASS[i][j] = MASS[i][j] + dt*resMass[i][j];
    XMOM[i][j] = XMOM[i][j] + dt*resXMom[i][j];
    YMOM[i][j] = YMOM[i][j] + dt*resYMom[i][j];
    ENER[i][j] = ENER[i][j] + dt*resEner[i][j];
  }
}


}
