#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argn, char** args){

			printf("Welcome to CFD Simulator for Arbitray Geometries:");
			printf("Worksheet 2, Group D: Moritz (Travis), Moe, Nia \n");
			printf("Please make a selection the scenerio you would like to simulate by a number between 1-7 \n");
			printf("P1. Karman Vortex Street \n");
			printf("P2. Flow over a Step \n");
			printf("P3. Natural Convection low reynolds number \n");
			printf("P4. Natural Convection High reynolds number\n" );
			printf("P5. Fluid Trap \n");
			printf("P6. Fluid Trap Reversed\n");
			printf("P7. Rayleigh-Benard Convection \n");
			int choice;
			char* geometry = (char*)(malloc(sizeof(char)*20));
			char* problem = (char*)(malloc(sizeof(char)*20));
			scanf("%d",&choice);

			const char* filename = "0";
			switch(choice)
			{
			case 1:
			filename = "Data_Files/karman_vortex.dat";

			break;
			case 2:
			filename = "Data_Files/step_flow.dat";

			break;
			case 3:
			filename = "Data_Files/natural_convection.dat";
			break;

			case 4:
			filename = "Data_Files/natural_convection2.dat";
			break;

			case 5:
			filename = "Data_Files/fluid_trap.dat";
			break;

			case 6:
			filename = "Data_Files/fluid_trap_reversed.dat";

			break;
			case 7:
			filename = "Data_Files/rb_convection.dat";
			break;

}


    double Re;
    double UI;
		double VI;
    double PI;
    double GX;
    double GY;
    double t_end;
    double xlength;
    double ylength;
    double dt;
    double dx;
    double dy;
    int    imax;
    int    jmax;
    double alpha;
    double omg;
    double tau;
    int    itermax;

    double eps;
    double dt_value;
    double Pr;
    double TI;
    double T_h;
    double T_c;
    double beta;

    //Read and assign the parameter values from file
    read_parameters(filename, &imax, &jmax, &xlength, &ylength,
			&dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax,
			&GX, &GY, &Re, &Pr, &UI, &VI, &PI, &TI, &T_h, &T_c, &beta, &dx, &dy, problem, geometry);


		int includeTemp = 1;
		if((choice == 1) || (choice == 2)){
			includeTemp=0;
		}


    double **P = matrix(0, imax-1, 0, jmax-1);
    double **U = matrix(0, imax-1, 0, jmax-1);
    double **V = matrix(0, imax-1, 0, jmax-1);
    double **F = matrix(0, imax-1, 0, jmax-1);
    double **G = matrix(0, imax-1, 0, jmax-1);
    double **RS = matrix(0, imax-1, 0, jmax-1);
    int **flag = imatrix(0, imax-1, 0, jmax-1);
    double **T;
    double **T1;

	if(includeTemp){
		T = matrix(0, imax-1, 0, jmax-1);
		T1= matrix(0, imax-1, 0, jmax-1);
		}


    //Initilize flags
		char g[100];
		 strcpy(g,"PGM_Files/");
		 strcat(g,geometry);
		 printf("%s\n", g);
    init_flag(problem, g, imax, jmax, flag);

    //Initialize the U, V and P
    if(includeTemp)
	{
		init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag);
	}
	else
	{
		init_uvp(UI, VI, PI, imax, jmax, U, V, P, flag);
	}

	// create a solution folder.



	 char temp[100];

	 strcpy (temp,"FINAL_");

	 strcat (temp,problem);

	 printf("%s\n",problem );

  double t=0; int n=0; int n1=0;

	while (t < t_end) {

		calculate_dt(Re,tau,&dt,dx,dy,imax,jmax, U, V, Pr, includeTemp);
   		//printf("t = %f ,dt = %f, ",t,dt);

		boundaryvalues(imax, jmax, U, V, flag,choice);

		if(includeTemp){
			calculate_temp(T, T1, Pr, Re, imax, jmax, dx, dy, dt, alpha, U, V, flag, TI, T_h, T_c, choice);
		}

    	spec_boundary_val(imax, jmax, U, V, flag);

    	calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G,flag, beta, T, includeTemp);

    	calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);

		int it = 0;
		double res = 10.0;

		struct stat st = {0};
		char solutionFolder[80];
		sprintf( solutionFolder,"Solution_%s",problem);
		if (stat(solutionFolder, &st) == -1) {
					mkdir(solutionFolder, 0700);
		}

		char solutionDirectory[80];
		sprintf( solutionDirectory,"Solution_%s/sol", problem);

    	do {
    		sor(omg,dx,dy,imax,jmax,P,RS,&res,flag);
				++it;

    	} while(it<itermax && res>eps);

			//printf("SOR itertions = %d ,residual = %f \n", it-1, res);
		  calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);

		if(includeTemp){
			remove_obstacles_temp(U, V, P, T, flag, imax, jmax);
  		}
		else{
			remove_obstacles(U, V, P, flag, imax, jmax);
			}

		if((t >= n1*dt_value)&&(t!=0.0)){
   			write_vtkFile(solutionDirectory,n ,xlength ,ylength ,imax-2 ,jmax-2 , dx ,dy ,U ,V ,P,T,includeTemp);

			printf("writing results for %f seconds \n",n1*dt_value);
    		++n1;
    		continue;
  		}

	//	if(n%32 == 0 )write_vtkFile(solutionDirectory, n, xlength, ylength, imax, jmax, dx, dy, U, V, P,T,includeTemp);
	//	if(n%250 == 0 ) printf("\nCurrent Time: %03f",t);
    	t =t+dt;
    	++n;
    }
			write_vtkFile(temp,n ,xlength ,ylength ,imax-2 ,jmax-2 , dx ,dy ,U ,V ,P,T,includeTemp);

    //Free memory
    free_matrix( P, 0, imax-1, 0, jmax-1);
    free_matrix( U, 0, imax-1, 0, jmax-1);
    free_matrix( V, 0, imax-1, 0, jmax-1);
    free_matrix( F, 0, imax-1, 0, jmax-1);
    free_matrix( G, 0, imax-1, 0, jmax-1);
    free_matrix(RS, 0, imax-1, 0, jmax-1);
    free_imatrix(flag, 0, imax-1, 0, jmax-1);
	if(includeTemp) { free_matrix(T, 0, imax-1, 0, jmax-1);
			   free_matrix(T1, 0, imax-1, 0, jmax-1); }
	free(geometry);
	free(problem);

	printf("End of Simulation.\n");
  return -1;

}
