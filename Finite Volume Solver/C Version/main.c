#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "solver.h"
#include <stdio.h>
#include <stdlib.h>


// void writeMatrix(int** m,int imax,int jmax,char* n)
// {
//    FILE * fp;
//    int i,j;
//
//    fp = fopen (n,"w");
//
//    for(j = jmax+1; j >= 0;j--){
//    for(i = 0; i <= imax+1;i++){
//        fprintf (fp, "%03d ",m[i][j]);
//    }
//    fprintf(fp,"\n");
//    }
//
//    /* close the file*/
//    fclose (fp);
//
// }

int main(int argn, char** args){
  printf("########### Welcome to the Finite Volume Solver for Computational Fluid Dynamics ##############\n");
  printf("Please make a selection for which scenerio you would like to simulate:\n");
  printf("1. Driven Cavity\n" );
  printf("2. Supersonic Compressible Flow\n");
  printf("3. Supersonic Compressible Flow Through Nozzle\n");
  printf("4. Shallow Water Equation\n");
  printf("5. Breaking Dam\n");
  int choice;

  scanf("%d", &choice);

  switch (choice){
    case 1:
      printf("Simulation starting for: Driven Cavity\n");
      cavity_fv();
      break;
    case 2:
      printf("Simulation starting for: Supersonic Compressible Flow\n");
      printf("You will be asked to choice your geomerty as well as a desired mach number\n" );
      supersonic_compressible();
      break;
    case 3:
      printf("Simulation of Supersonic Flow Through a Nozzle\n");
      supersonic_nozzle();
    case 4:
      printf("Simulation starting for: Shallow Water\n");
      shallow_water();
      break;
    case 5:
      printf("Dam Break Simulation\n");
      breaking_dam();
      break;
    default:
      printf("Choice not valid\n");
      break;
  }

  return -1;
}
