#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
#include <stdio.h>

int main(int argc, char* argv[])
{
    const char* filename= "/Users/moritzhof/Documents/TUM/Semester2/ComputationalFluidDynamicsLab/Worksheet 1/cavity100.dat";
    //const char* destination = "/Users/moritzhof/Documents/TUM/Semester 2/ComputationalFluidDynamicsLab/Worksheet 1/results";
	//initialize variables
	double t = 0; /*time start*/
	int it, n = 0; /*iteration and time step counter*/
	double res; /*residual for SOR*/
		/*arrays*/
	double **U, **V, **P;
	double **RS, **F, **G;
		/*those to be read in from the input file*/
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value;
	int  imax, jmax, itermax;


	//read the parameters
    read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax,
        &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value);

	//allocate memory
	U = matrix(0, imax+1, 0, jmax+1);
	V = matrix(0, imax+1, 0, jmax+1);
	P = matrix(0, imax+1, 0, jmax+1);
	RS = matrix(1, imax, 1, jmax);
	F = matrix(0, imax, 1, jmax);
	G = matrix(1, imax, 0, jmax);


	init_uvp(UI, VI, PI, imax, jmax, U, V, P);
  printf("dt=%f\n", dt);
	//going through all time steps
	while(t < t_end){
		//adaptive time stepping
		//calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);

		//setting bound.values
		boundaryvalues(imax, jmax, U, V);

		//computing F, G and right hand side of pressue eq.
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		//iteration counter
		it = 0;

		do{
			//perform SOR iteration, at same time set bound.values for P and new residual value
			sor(omg, dx, dy, imax, jmax, P, RS, &res);
			it++;
		}while(it<itermax && res>eps);
		//calculate U and V of this time step
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
    //printf("%f\t", U[(1/2)*imax][7/8*jmax]);
    // if (n%4==0){
		// 	write_vtkFile("solutions", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
		// }
		//indent time and number of time steps
		t += dt;
    n++;
		//output of pics for animation
	}
    printf("dt=%f\n", dt);
    //printf("imax = %d, jmax = %d \n", (1/2)*imax, (7*8)*jmax);
    printf("U=%f\n", U[imax/2][(7/8*jmax)]);
	//output of U, V, P at the end for visualization
  write_vtkFile("solutions_default_settings", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

	//free memory
	free_matrix(U, 0, imax+1, 0, jmax+1);
	free_matrix(V, 0, imax+1, 0, jmax+1);
	free_matrix(P, 0, imax+1, 0, jmax+1);
	free_matrix(RS, 1, imax, 1, jmax);
	free_matrix(F, 0, imax, 1, jmax);
	free_matrix(G, 1, imax, 0, jmax);

	return -1;
}
