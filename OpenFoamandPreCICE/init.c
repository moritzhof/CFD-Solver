#include "helper.h"
#include "init.h"
#include <stdio.h>

void read_parameters( const char *szFileName,
				char* problem,
				char* geometry,
	  		char* precice_config,
				char* participant_name,
				char* mesh_name,
				char* read_data_name,
				char* write_data_name,
		    int  *imax,
        int  *jmax,
		    double *xlength,
        double *ylength,
				double *x_origin,
				double *y_origin,
		    double *dt,
		    double *t_end,
		    double *tau,
		    double *dt_value,
		    double *eps,
		    double *omg,
		    double *alpha,
        int  *itermax,
		    double *GX,
        double *GY,
		    double *Re,
        double *Pr,
		    double *UI,
        double *VI,
        double *PI,
        double *TI,
		    double *beta,
		    double *dx,
        double *dy

)
{
   printf("PROGRESS: Reading .dat file... \n");
   READ_STRING( szFileName, problem );
   READ_STRING( szFileName, geometry );

	 READ_STRING(szFileName, precice_config);
	 READ_STRING(szFileName, participant_name);
	 READ_STRING(szFileName, mesh_name);
	 READ_STRING(szFileName, read_data_name);
	 READ_STRING(szFileName, write_data_name);

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );
	 READ_DOUBLE( szFileName, *x_origin );
	 READ_DOUBLE( szFileName, *y_origin);

   READ_DOUBLE( szFileName, *dt    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *alpha );
   READ_INT   ( szFileName, *itermax );

   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *Re );
   READ_DOUBLE( szFileName, *Pr );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *TI );
   READ_DOUBLE( szFileName, *beta );


   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   printf("Parameters: .dat read in \n");

}


void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax,
		 					double** U, double** V, double** P, double** T, int** flag){

	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if(flag[i][j]&(1<<0)){
					U[i][j] = UI;
					V[i][j] = VI;
					P[i][j] = PI;
					T[i][j] = TI;
			}
		}
	}
	printf("INIT_UVPT: Done\n");
}

int fluid(int pic){
	if((pic == 2)||(pic == 3)||(pic == 6)) {return 1;}
		else {return 0;}
}



void init_flag(char* problem, char* geometry, int imax, int jmax, int **flag, int *num_coupling_cells){

	printf("PROGRESS: Setting flags... \n");

	int **pic = imatrix(0,imax-1,0,jmax-1);
	pic = read_pgm(geometry);

	//initialiation: Coupling Cells:
	*num_coupling_cells = 0;

	for (int i=0; i<imax; i++){
		for (int j=0; j<jmax; j++){

		flag[i][j] = 0;


		switch(pic[i][j]){
			case 0:
			flag[i][j] = 1<<1;
			break;

			case 1:
			flag[i][j] = 1<<2;
			break;

			case 2:
			flag[i][j] = 1<<3;
			break;

			case 3:
			flag[i][j] = 1<<4;
			break;

			case 4: //coupling | no-slip
			flag[i][j] = 1<<9 | 1<<1;
			(*num_coupling_cells)++;
			break;

			case 6:
			 	flag[i][j] = 1<<0;
			 break;
		}

			if(!fluid(pic[i][j])){
				if(i<imax-1 && pic[i+1][j]==6){
					flag[i][j] |= 1<<8;
				}

				if(i>0 && pic[i-1][j]==6){
					flag[i][j] |= 1<<7;
				}
				if(j<jmax-1 && pic[i][j+1]==6){
					flag[i][j] |= 1<<5;
				}
				if(j>0 && pic[i][j-1]==6){
					flag[i][j] |= 1<<6;
				}
			}
		}
	}
	printf("Number of coupled cells: %d\n", *num_coupling_cells );
	free_imatrix(pic, 0,imax-1,0,jmax-1);
	printf("INIT FLAG: DONE \n");


}
