#include"helper.h"
#include"init.h"
#include"uvp.h"
#include"visual.h"

#include<iostream>
#include<iterator>
#include<vector>
#include<math.h>
#include<numeric>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include<utility>

//##############################################################################
//################# C++14 Version and above ####################################
//##############################################################################

int main(int argc, char const *argv[]){

//########################## Initialization ###################################

  int imax = 500, jmax = 200;
  int dimension = imax*jmax;
  double rho = 1.225, p = 101325.0;
  double mach = 2.5, k = 1.4, ux = mach*std::sqrt(k*p/rho);
  double uy = 0;
  double enerConst = 0.5*rho*(ux*ux+uy*uy)+(p/(k-1));
  double xlength = 5, ylength = 2;

  double cfl = 0.6;
  double time = 0, t_end = 0.3, dt = 0.05;

  double dx = xlength/imax, dy = ylength/jmax;

  std::cout<<enerConst<<std::endl;

  std::cout<<  0.5*rho*(ux*ux+uy*uy)+(p/(k-1)) << std::endl;
//########################## Matrices && Initialization ########################

  //Initialization
  Matrix MASS(imax+2, std::vector<double>(jmax+2));
  Matrix XMOM(imax+2, std::vector<double>(jmax+2));
  Matrix YMOM(imax+2, std::vector<double>(jmax+2));
  Matrix ENER(imax+2, std::vector<double>(jmax+2));


  Matrix resMass(imax+2, std::vector<double>(jmax+2));
  Matrix resXMom(imax+2, std::vector<double>(jmax+2));
  Matrix resYMom(imax+2, std::vector<double>(jmax+2));
  Matrix resEner(imax+2, std::vector<double>(jmax+2));

  intMatrix FLAG(imax+2, std::vector<int>(jmax+2));
  intMatrix PIC(imax+2, std::vector<int>(jmax+2));

  //Initialization to a non-zero value;
  init_matrix(MASS, rho);
  init_matrix(XMOM, rho*ux);
  init_matrix(YMOM, rho*uy);
  init_matrix(ENER, enerConst);

  //Initialization of .pmg
  const char* geometry = "rocket.pgm";
  //intMatrix Pic = read_pgm(geometry);

  //std::cout<< Pic.size()<< " "<< Pic[0].size()<<std::endl;

  init_FLAG(geometry, FLAG, PIC);
  std::cout<<" FLAG Done" <<std::endl;

//   for(int i = 0; i<imax+1; ++i){
//     for(int j = 0; j<jmax+1; ++j){
//     printf("%f ", ENER[i][j]);
//   }
//   printf("\n");
// }

  //printMatrix(FLAG);
  int init_value = 0;
  int n = 0;
  while(time < t_end){
      init_matrix(resMass, init_value);
      init_matrix(resXMom, init_value);
      init_matrix(resYMom, init_value);
      init_matrix(resEner, init_value);


  if(n%50 == 0) {
  write_vtkFile("Solution",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);}
  xLocalLaxFriedrichs(&dt,imax,jmax,rho,ux,uy,p,k,dx,cfl,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner, FLAG);
  yLocalLaxFriedrichs(&dt,imax,jmax,rho,ux,uy,p,k,dy,cfl,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner, FLAG);
  timeStep(dt,imax,jmax,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner);
  //printMatrix(MASS);
  time += dt;
  ++n;
}
  std::cout<<"Simulation Complete"<<std::endl;
  //printMatrix(YMOM);
return 0;
}
