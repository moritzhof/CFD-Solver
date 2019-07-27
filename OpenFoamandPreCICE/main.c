#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "precice_adapter.h"
#include "precice/SolverInterfaceC.h"


int main(int argc, char const *argv[]){

			printf("Welcome to CFD Simulator with PreCICE and OpenFoam:");
			printf("Please make a selection the scenerio you would like to simulate by a number between 1-7 \n");
			printf("1. Heated Plate: Forced Convection \n");
			printf("2. Conducting Wall: Natural Convection\n");
			printf("3. Heat Exchanger: Fluid 1 \n");
			printf("4. Heat Exchanger: Fluid 2\n");

			int choice;
			char* geometry = (char*)(malloc(sizeof(char)*256*4));
			char* problem = (char*)(malloc(sizeof(char)*256*4));


			//These to be added to read_parameters:
			//precice_config.xml
			char* precice_config =(char*)(malloc(sizeof(char)*500));
			//Fluid
			char* participant_name =(char*)(malloc(sizeof(char)*500));
			//Fluid Mesh
			char* mesh_name = (char*)(malloc(sizeof(char)*500));
			//Heat Hlux
			char* read_data_name = (char*)(malloc(sizeof(char)*500));
			//Temperature
			char* write_data_name = (char*)(malloc(sizeof(char)*500));


			scanf("%d",&choice);

			const char* filename = "0";
			switch(choice){

			case 1:
				filename = "configs/heated-plate.dat";
			break;

			case 2:
				filename = "configs/convection.dat";
			break;

			case 3:
				filename = "configs/F1-heat-exchange.dat";
			break;

			case 4:
				filename = "configs/F2-heat-exchange.dat";
			break;

			}

			const char* coric = precicec_actionReadIterationCheckpoint();
      const char* cowic = precicec_actionWriteIterationCheckpoint();



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
		double time_cp;
		double x_origin, y_origin;

		read_parameters(filename, problem, geometry, precice_config, participant_name, mesh_name, read_data_name, write_data_name, &imax, &jmax, &xlength, &ylength,
										&x_origin, &y_origin, &dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax, &GX, &GY, &Re, &Pr, &UI, &VI, &PI, &TI, &beta, &dx,
										&dy);

    double **P = matrix(0, imax, 0, jmax);
    double **U = matrix(0, imax, 0, jmax);
    double **V = matrix(0, imax, 0, jmax);
    double **F = matrix(0, imax, 0, jmax);
    double **G = matrix(0, imax, 0, jmax);
    double **RS = matrix(0, imax, 0, jmax);
    int **flag = imatrix(0, imax, 0, jmax);
    double **T = matrix(0, imax, 0, jmax);
		double ** T_cp = matrix(0, imax, 0, jmax);
    double **T1 = matrix(0, imax, 0, jmax);
		double **U_cp = matrix(0, imax, 0, jmax);
		double **V_cp = matrix(0, imax, 0, jmax);



		//Initilize preCICE:
		precicec_createSolverInterface(participant_name, precice_config, 0,1);
		int dim = precicec_getDimensions();

	   int meshID = precicec_getMeshID(mesh_name);
		 int num_coupling_cells = 0;


		 //Initilize flags
     init_flag(problem, geometry, imax, jmax, flag, &num_coupling_cells);

		 init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag);

		 int* vertexIDs = precice_set_interface_vertices(imax, jmax, dx, dy, x_origin, y_origin, num_coupling_cells, meshID, flag);
		 printf("VextexIDs: done\n");
		//define Dirichlet part of the couple to be written by this solver
		int temperatureID = precicec_getDataID(write_data_name, meshID);

		double* temperatureCoupled = (double*)malloc(sizeof(double)*num_coupling_cells);

		//define Neuman part of the coupling read by this solver.
		int heatFluxID = precicec_getDataID(read_data_name, meshID);
		double* heatFluxCoupled = (double*)malloc(sizeof(double)*num_coupling_cells);
		printf("Temp and heat flux: done\n");

		//call precicec_initialize()
		double precice_dt = precicec_initialize();


		precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled, vertexIDs, temperatureID, T, flag);
		//openFoam
		precicec_initialize_data();
		precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatFluxCoupled);

	// create a solution folder.
	struct stat st = {0};
	char solutionFolder[80];
	sprintf( solutionFolder,"%s",problem);
	if (stat(solutionFolder, &st) == -1) {
    		mkdir(solutionFolder, 0700);
	}

	char solutionDirectory[80];
	sprintf(solutionDirectory,"%s/sol", problem);


	int n = 0, n1 =0;
	double t = 0; double res = 10.0;
	while (precicec_isCouplingOngoing()){

			if(precicec_isActionRequired(cowic)){
				write_checkpoint(t, U, V, T, &time_cp, U_cp, V_cp, T_cp, imax, jmax);
				precicec_fulfilledAction(cowic);
			}

			calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V, Pr);
			dt = fmin(dt, precice_dt);
			boundaryvalues(imax, jmax, U, V, flag, UI);
			set_coupling_boundary(imax,jmax, dx, dy, heatFluxCoupled, T,  flag);
			calculate_temp(T, T1, Pr, Re, imax, jmax, dx, dy, dt, alpha, U , V, flag, TI, UI);
			spec_boundary_val(imax, jmax, U, V, flag, choice);
			calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, flag, beta, T, UI);
			calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, flag);

			int it = 0;
			double res = 10.0;
			do{
				sor(omg, dx, dy, imax, jmax, P, RS, &res, flag, UI);
				++it;
			}while (it<itermax && res>eps);

			calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, flag, UI);

			precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled, vertexIDs, temperatureID, T, flag);
			precice_dt = precicec_advance(dt);
			precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatFluxCoupled);


			if ((t >= n1*dt_value)&&(t!=0.0)){
				write_vtkFile(solutionDirectory ,n1 ,xlength ,ylength, x_origin, y_origin,imax-1 ,jmax-1 ,dx ,dy ,U ,V ,P, T);
				++n1;
				continue;
			}
			remove_obstacles_temp(U, V, P, T, flag, imax, jmax);
			if(precicec_isActionRequired(coric)){
				restore_checkpoint(&t, U, V, T, time_cp, U_cp, V_cp, T_cp, imax, jmax);
				precicec_fulfilledAction(coric);
				}
				else{
					t=t+dt;
					++n;
			}


    }
		precicec_finalize();

    //Free memory
    free_matrix( P, 0, imax, 0, jmax);
    free_matrix( U, 0, imax, 0, jmax);
    free_matrix( V, 0, imax, 0, jmax);
    free_matrix( F, 0, imax, 0, jmax);
    free_matrix( G, 0, imax, 0, jmax);
    free_matrix(RS, 0, imax, 0, jmax);
    free_imatrix(flag, 0, imax, 0, jmax);
	  free_matrix(T, 0, imax, 0, jmax);
		free_matrix(T1, 0, imax, 0, jmax);
		free_matrix(U_cp, 0, imax, 0, jmax);
    free_matrix(V_cp, 0, imax, 0, jmax);
	  free(geometry);
	  free(problem);
		free(precice_config);
		free(participant_name);
		free(mesh_name);
		free(read_data_name);
		free(write_data_name);
		free(vertexIDs);
		free(temperatureCoupled);
		free(heatFluxCoupled);

	printf("End of Simulation.\n");
  return -1;

}
