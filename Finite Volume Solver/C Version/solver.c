#include "solver.h"
#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "boundary_val.h"
#include <stdio.h>

void supersonic_nozzle(){

  int choice;
  char* geometry;
  printf("Please Choose the Type of Geometry\n");
  printf("1. Air Spike Nozzle\n");
  printf("2. Cone Nozzle\n");
  printf("3. Bell Nozzle\n");

  scanf("%i", &choice);
  switch (choice) {
    case 1:
      geometry = "nozzle_airspike.pgm";
      break;
    case 2:
      geometry = "nozzle_cone.pgm";
      break;
    case 3:
      geometry = "b_n.pgm";
      break;
    default:
      printf("Not a Valid Choice\n");
      break;
  }
  int imax = 2000;
  int jmax = 400;
  double xlength =  20.0;
  double ylength =  4.0;
  int **flag;

  flag = imatrix(0, imax+1, 0, jmax+1);
  init_flag(geometry,imax,jmax,flag);
  // writeMatrix(flag,imax,jmax,"1.txt");


//GET THESE VALUES FROM read_parameters
  double rho    = 1.225;
  double p      = 101325.0;
  double mach   = 0.0;
  double k      = 1.4;
  double ux     = mach*sqrt(k*p/rho);
  double uy     = 0.0;





  double cfl = 0.5;
  double t_end = 1.0;
  double time = 0.0;
  double dt = 0.05;

  double dx = xlength/imax;
  double dy = ylength/jmax;
  printf("%f %f\n",dx,dy );

  //Properties that are conserved (hyperbolic terms)
  double** MASS = matrix(0,imax+1,0,jmax+1);
  double** XMOM = matrix(0,imax+1,0,jmax+1);
  double** YMOM = matrix(0,imax+1,0,jmax+1);
  double** ENER = matrix(0,imax+1,0,jmax+1);

  //Residual Matrices (DUdt) for each terms
  double** resMass = matrix(0,imax+1,0,jmax+1);
  double** resXMom = matrix(0,imax+1,0,jmax+1);
  double** resYMom = matrix(0,imax+1,0,jmax+1);
  double** resEner = matrix(0,imax+1,0,jmax+1);

  init_ssc(rho,0,0, k, p,imax,jmax,MASS,XMOM,YMOM,ENER);

  for(int i=0;i<125;++i){
    for(int j=0;j<400;++j){
      if(flag[i][j] & FLUID){
        MASS[i][j] = 10.0*rho;
        ENER[i][j] = 0.5*rho*(ux*uy)+100.0*p/(k-1);
      }
    }
  }

int n = 0;


while(time < t_end){
  init_matrix(resMass,0,imax+1,0,jmax+1,0);
  init_matrix(resXMom,0,imax+1,0,jmax+1,0);
  init_matrix(resYMom,0,imax+1,0,jmax+1,0);
  init_matrix(resEner,0,imax+1,0,jmax+1,0);


  if(n%500 == 0){
    if(choice == 1){
    write_vtkFile_ssc("air_spike",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
    }
    if(choice == 2){
      write_vtkFile_ssc("cone_nozzle",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
    }
    if(choice == 3){
      write_vtkFile_ssc("bell_nozzle",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
    }
  }
    boundaryvalues_nozzle(imax,jmax,rho,ux,uy,p,k,MASS,XMOM,YMOM,ENER,flag);


  x_llf_ssc(&dt,imax,jmax,rho,ux,uy,p,k,dx,cfl,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner,flag);
  y_llf_ssc(&dt,imax,jmax,rho,ux,uy,p,k,dy,cfl,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner,flag);
// if(n/10 == 500)break;
  time_step_ssc(dt,imax,jmax,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner);

  time  += dt;


  n++;



}


  free_imatrix(flag,0, imax+1, 0, jmax+1);

  free_matrix(MASS,0, imax+1, 0, jmax+1);
  free_matrix(YMOM,0, imax+1, 0, jmax+1);
  free_matrix(XMOM,0, imax+1, 0, jmax+1);
  free_matrix(ENER,0, imax+1, 0, jmax+1);

  free_matrix(resMass,0, imax+1, 0, jmax+1);
  free_matrix(resXMom,0, imax+1, 0, jmax+1);
  free_matrix(resYMom,0, imax+1, 0, jmax+1);
  free_matrix(resEner,0, imax+1, 0, jmax+1);



}

//##############################################################################################
//##############################################################################################
void supersonic_compressible()
{


  char* geometry = "";
  printf(" Please select the shape:\n");
  printf("1. Rocket\n");
  printf("2. Cylinder\n");
  int choice;
  scanf("%i", &choice);
  switch (choice){
    case 1:
      geometry = "rocket.pgm";
      break;
  case 2:
      geometry = "rocket2.pgm";
      break;
  default:
      printf("Not a valid \n");
      break;
  }

  int imax = 500;
  int jmax = 200;

  int **flag;
  flag = imatrix(0, imax+1, 0, jmax+1);
  init_flag(geometry,imax,jmax,flag);



//GET THESE VALUES FROM read_parameters
  double rho    = 1.225;
  double p      = 101325.0;
  double mach;//   = 2.5;
  printf("Please choice a mach number greater than 1\n");
  scanf("%lf", &mach);
  double k      = 1.4;
  double ux     = mach*sqrt(k*p/rho);
  double uy     = 0;



  double xlength =  5;
  double ylength =  2;

  double cfl = 0.6;
  double t_end = 0.3;
  double time = 0;
  double dt = 0.05;
  // double dt_value = 0.01;
  // double writeTime = 0;

  double dx = xlength/imax;
  double dy = ylength/jmax;
  printf("%f %f\n",dx,dy );

  //Properties that are conserved (hyperbolic terms)
  double** MASS = matrix(0,imax+1,0,jmax+1);
  double** XMOM = matrix(0,imax+1,0,jmax+1);
  double** YMOM = matrix(0,imax+1,0,jmax+1);
  double** ENER = matrix(0,imax+1,0,jmax+1);

  //Residual Matrices (DUdt) for each terms
  double** resMass = matrix(0,imax+1,0,jmax+1);
  double** resXMom = matrix(0,imax+1,0,jmax+1);
  double** resYMom = matrix(0,imax+1,0,jmax+1);
  double** resEner = matrix(0,imax+1,0,jmax+1);

  init_ssc(rho,ux,uy, k, p,imax,jmax,MASS,XMOM,YMOM,ENER);
int n = 0;

write_vtkFile_ssc("Solution",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
while(time < t_end){
  init_matrix(resMass,0,imax+1,0,jmax+1,0);
  init_matrix(resXMom,0,imax+1,0,jmax+1,0);
  init_matrix(resYMom,0,imax+1,0,jmax+1,0);
  init_matrix(resEner,0,imax+1,0,jmax+1,0);

  // if(time - writeTime >= dt_value){
  if(n%10 == 0){
    printf("Writing into file at time = %f dt = %f\n",time,dt );
    if(choice == 1)
      write_vtkFile_ssc("Rocket_Solution",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
    else if(choice == 2)
      write_vtkFile_ssc("Cylinder_Solution",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
    //write_vtkFile_ssc("Solution",n,xlength,ylength,imax,jmax,dx,dy,k,MASS,XMOM,YMOM,ENER);
    //writeTime = time;
}

  dt = 0.5;
  x_llf_ssc(&dt,imax,jmax,rho,ux,uy,p,k,dx,cfl,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner,flag);
  y_llf_ssc(&dt,imax,jmax,rho,ux,uy,p,k,dy,cfl,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner,flag);

  time_step_ssc(dt,imax,jmax,MASS,XMOM,YMOM,ENER,resMass,resXMom,resYMom,resEner);

  time  += dt;


  n++;



}

  free_imatrix(flag,0, imax+1, 0, jmax+1);

  free_matrix(MASS,0, imax+1, 0, jmax+1);
  free_matrix(YMOM,0, imax+1, 0, jmax+1);
  free_matrix(XMOM,0, imax+1, 0, jmax+1);
  free_matrix(ENER,0, imax+1, 0, jmax+1);

  free_matrix(resMass,0, imax+1, 0, jmax+1);
  free_matrix(resXMom,0, imax+1, 0, jmax+1);
  free_matrix(resYMom,0, imax+1, 0, jmax+1);
  free_matrix(resEner,0, imax+1, 0, jmax+1);

}

//###################################################################################
void cavity_fv()
{
  int imax = 50;
  int jmax = 50;
  double xlength =  1;
  double ylength =  1;

  int **flag;

  flag = imatrix(0, imax+1, 0, jmax+1);

init_flag("cavity.pgm",imax,jmax,flag);

//GET THESE VALUES FROM read_parameters


  double t_end = 10;
  double time = 0;
  double dt = 0.5;



  double dx = xlength/imax;
  double dy = ylength/jmax;
  printf("%f %f\n",dx,dy );


  double** U = matrix(0,imax+1,0,jmax+1);
  double** V = matrix(0,imax+1,0,jmax+1);
  double** P    = matrix(0,imax+1,0,jmax+1);
  double** RS   = matrix(1,imax  ,1,jmax  );



  double** F = matrix(0,imax+1,0,jmax+1);
  double** G = matrix(0,imax+1,0,jmax+1);

  double** resP = matrix(0,imax+1,0,jmax+1);


  init_uvp(0,0,0,imax,jmax,U,V,P);
  int n = 0;


double eps = 0.01;
double omg = 1.7;
double itmax = 1000;
double Re =100;
double tau = 0.5;
while(time < t_end){

  boundaryvalues(imax, jmax, U, V);
  calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);

if (n%50 == 0) { printf("time  = %f dt  = %f\n", time,dt);
  write_vtkFile_cavity("cavity",n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
}


  calculate_fg_cds_fv(dt,imax,jmax,dx,dy,U,V,F,G,flag);


  caluculate_p_fv(imax,jmax,itmax,dt,dx,dy,omg,eps,F,G,P,RS,resP);
  calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);

  printf("time  = %f dt  = %f\n", time,dt);
//
  time  += dt;
  n++;
//
//
}

write_vtkFile_cavity("cavity",n,xlength,ylength,imax,jmax,dx,dy,U,V,P);

  free_imatrix(flag,0, imax+1, 0, jmax+1);


  free_matrix(V,0, imax+1, 0, jmax+1);
  free_matrix(U,0, imax+1, 0, jmax+1);

  free_matrix(resP,0,imax+1,0,jmax+1);

  free_matrix(RS,1,imax,1,jmax);

  free_matrix(F,0, imax+1, 0, jmax+1);
  free_matrix(G,0, imax+1, 0, jmax+1);
}


//##################################################################################
void shallow_water(){
  int imax = 100;
  int jmax = 100;

//GET THESE VALUES FROM read_parameters
  double ux = 0;
  double uy = 0;
  double g  = 9.81;


  int **flag = imatrix(0, imax+1, 0, jmax+1);
  init_flag_sw("droplet.pgm",imax,jmax,flag);

  double xlength =  1;
  double ylength =  1;

  double cfl = 0.5;
  double t_end = 2;
  double time = 0;
  double dt = 0.005;

  double dx = xlength/imax;
  double dy = ylength/jmax;
  printf("%f %f\n",dx,dy );

  //Properties that are conserved (hyperbolic terms)
  double** MASS = matrix(0,imax+1,0,jmax+1);
  double** XMOM = matrix(0,imax+1,0,jmax+1);
  double** YMOM = matrix(0,imax+1,0,jmax+1);

  //Residual Matrices (DUdt) for each terms
  double** resMass = matrix(0,imax+1,0,jmax+1);
  double** resXMom = matrix(0,imax+1,0,jmax+1);
  double** resYMom = matrix(0,imax+1,0,jmax+1);


  int base;
  double height;
  printf("Please enter an Base and a Height. Note: Base must be greater than 0 and of type int, and max height 10\n");
  printf("Please enter the base: \n");
  scanf("%i", &base);
  printf("Please enter the height: \n");
  scanf("%lf", &height);


  init_sw(base,ux,uy,imax,jmax,MASS,XMOM,YMOM);

  //Inintial hieght of wter
    for(int i = 45; i<=55; i++){
      for(int j = 45; j<=55; j++){
        // if(flag[i][j] & OUTFLOW){
          MASS[i][j] = height;
          XMOM[i][j] = height*ux;
          YMOM[i][j] = height*uy;
        }
      }

int n = 0;

write_vtkFile_sw("shallow_water",n,xlength,ylength,imax,jmax,dx,dy,MASS);
while(time < t_end){
  init_matrix(resMass,0,imax+1,0,jmax+1,0);
  init_matrix(resXMom,0,imax+1,0,jmax+1,0);
  init_matrix(resYMom,0,imax+1,0,jmax+1,0);

if(n%10 == 0){
    printf("Writing into file at time = %f\n",time );
    write_vtkFile_sw("shallow_water",n,xlength,ylength,imax,jmax,dx,dy,MASS);
}
  boundaryvalues_sw(imax,jmax,MASS,XMOM,YMOM,flag);
  x_llf_sw(&dt,imax,jmax,g,dx,cfl,MASS,XMOM,YMOM,resMass,resXMom,resYMom,flag);
  y_llf_sw(&dt,imax,jmax,g,dy,cfl,MASS,XMOM,YMOM,resMass,resXMom,resYMom,flag);

  time_step_sw(dt,imax,jmax,MASS,XMOM,YMOM,resMass,resXMom,resYMom,flag);

  time  += dt;


  n++;

}

  free_imatrix(flag,0, imax+1, 0, jmax+1);


  free_matrix(MASS,0, imax+1, 0, jmax+1);
  free_matrix(YMOM,0, imax+1, 0, jmax+1);
  free_matrix(XMOM,0, imax+1, 0, jmax+1);

  free_matrix(resMass,0, imax+1, 0, jmax+1);
  free_matrix(resXMom,0, imax+1, 0, jmax+1);
  free_matrix(resYMom,0, imax+1, 0, jmax+1);
}

void writeMatrix(int** m,int imax,int jmax,char* n)
{
   FILE * fp;
   int i,j;

   fp = fopen (n,"w");

   for(j = jmax+1; j >= 0;j--){
   for(i = 0; i <= imax+1;i++){
       fprintf (fp, "%03d ",m[i][j]);
   }
   fprintf(fp,"\n");
   }

   /* close the file*/
   fclose (fp);

}

//#############################################################################
void breaking_dam(){

  int imax = 100;
  int jmax = 300;

//GET THESE VALUES FROM read_parameters
  double ux = 0;
  double uy = 0;
  double g  = 9.81;


    int **flag = imatrix(0, imax+1, 0, jmax+1);
     init_flag_sw("dam.pgm",imax,jmax,flag);

  double xlength =  1;
  double ylength =  3.;

  double cfl = 0.5;
  double t_end = 2;
  double time = 0;
  double dt = 0.005;

  double dx = xlength/imax;
  double dy = ylength/jmax;
  printf("%f %f\n",dx,dy );

  //Properties that are conserved (hyperbolic terms)
  double** MASS = matrix(0,imax+1,0,jmax+1);
  double** XMOM = matrix(0,imax+1,0,jmax+1);
  double** YMOM = matrix(0,imax+1,0,jmax+1);

  //Residual Matrices (DUdt) for each terms
  double** resMass = matrix(0,imax+1,0,jmax+1);
  double** resXMom = matrix(0,imax+1,0,jmax+1);
  double** resYMom = matrix(0,imax+1,0,jmax+1);





  int base;
  double height;
  printf("Please enter an Base and a Height. Note: Base must be greater than 0 and of type int, and max height 10\n");
  printf("Please enter the base: \n");
  scanf("%i", &base);
  printf("Please enter the height: \n");
  scanf("%lf", &height);


  init_sw(base,ux,uy,imax,jmax,MASS,XMOM,YMOM);
  //Inintial hieght of wter
    for(int i = 0; i<=imax+1; i++){
      for(int j = 0; j<=jmax+1; j++){
        if(flag[i][j] & OUTFLOW){
          MASS[i][j] = height;
          XMOM[i][j] = height*ux;
          YMOM[i][j] = height*uy;
        }
      }
    }
int n = 0;


while(time < t_end){
  init_matrix(resMass,0,imax+1,0,jmax+1,0);
  init_matrix(resXMom,0,imax+1,0,jmax+1,0);
  init_matrix(resYMom,0,imax+1,0,jmax+1,0);

if(n%10 == 0){
    printf("Writing into file at time = %f\n",time );
write_vtkFile_sw("breaking_dam",n,xlength,ylength,imax,jmax,dx,dy,MASS);
}
boundaryvalues_sw(imax,jmax,MASS,XMOM,YMOM,flag);

  x_llf_sw(&dt,imax,jmax,g,dx,cfl,MASS,XMOM,YMOM,resMass,resXMom,resYMom,flag);
  y_llf_sw(&dt,imax,jmax,g,dy,cfl,MASS,XMOM,YMOM,resMass,resXMom,resYMom,flag);

  time_step_sw(dt,imax,jmax,MASS,XMOM,YMOM,resMass,resXMom,resYMom,flag);

  time  += dt;


  n++;


}

  free_imatrix(flag,0, imax+1, 0, jmax+1);


  free_matrix(MASS,0, imax+1, 0, jmax+1);
  free_matrix(YMOM,0, imax+1, 0, jmax+1);
  free_matrix(XMOM,0, imax+1, 0, jmax+1);

  free_matrix(resMass,0, imax+1, 0, jmax+1);
  free_matrix(resXMom,0, imax+1, 0, jmax+1);
  free_matrix(resYMom,0, imax+1, 0, jmax+1);

}
