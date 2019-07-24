#include "helper.h"
#include "init.h"
#include <stdio.h>

void read_parameters( const char *szFileName,    
		    int  *imax,
        int  *jmax,
		    double *xlength,
        double *ylength,
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
		    double *T_h,
		    double *T_c,
		    double *beta,
		    double *dx,
        double *dy,
        char *problem,
		    char *geometry

)
{
   printf("PROGRESS: Reading .dat file... \n");
   //READ_STRING( szFileName, *problem );
   //READ_STRING( szFileName, geometry );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

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
   READ_DOUBLE( szFileName, *T_h );
   READ_DOUBLE( szFileName, *T_c );
   READ_DOUBLE( szFileName, *beta );

	READ_STRING( szFileName, problem);
	READ_STRING( szFileName, geometry);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   printf("PROGRESS: .dat file read... \n \n");

}


void init_uvp(double UI, double VI, double PI, int imax, int jmax,
		 double** U, double** V, double** P, int** flag){

	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if(flag[i][j]&(1<<0)){

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
			}
		}
	}
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
}

int  fluid(int pic){
	if((pic == 2)||(pic == 3)||(pic == 4)) {return 1;}
		else {return 0;}
}



void init_flag(char* problem, char* geometry, int imax, int jmax, int **flag){

	int **pic = imatrix(0,imax-1,0,jmax-1);
	pic = read_pgm(geometry);

	for (int i=0; i<imax; i++){
		for (int j=0; j<jmax; j++){

		flag[i][j] = 0;


		switch(pic[i][j])
		{
			case 0: //fluid
			flag[i][j] = 1<<1;
			break;

			case 1: //no-slip
			flag[i][j] = 1<<2;
			break;

			case 2: //free-slip
			flag[i][j] = 1<<3;
			break;

			case 3: //outflow
			flag[i][j] = 1<<4;
			break;

			case 4: //inflow
			flag[i][j] = 1<<0;
			break;
		}

			if(!fluid(pic[i][j])){

				if(i<imax-1 && pic[i+1][j]==4)
				{
				flag[i][j] |= 1<<8;
				}
				if( i>0 && pic[i-1][j]==4)
				{
				flag[i][j] |= 1<<7;
				}
				if(j<jmax-1 && pic[i][j+1]==4)
				{
				flag[i][j] |= 1<<5;
				}
				if(j>0 && pic[i][j-1]==4)
				{
				flag[i][j] |= 1<<6;
				}
			}


		}

	}
	free_imatrix(pic, 0,imax-1,0,jmax-1);


}
