#include "helper.h"
#include "init.h"

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		    double *dt_value)           /* time for output */
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}

void init_uvp(
  double ux,
  double uy,
  double PI,
  int imax,
  int jmax,
  double **XMOM,
  double **YMOM,
  double ** P
)
{
  init_matrix(XMOM,0,imax+1,0,jmax+1,ux);
  init_matrix(YMOM,0,imax+1,0,jmax+1,uy);
  init_matrix(YMOM,0,imax+1,0,jmax+1,PI);
}

void init_ssc(
  double rho,
  double ux,
  double uy,
  double  k,
  double  p,
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER
)
{
  init_matrix(MASS,0,imax+1,0,jmax+1,rho);
  init_matrix(XMOM,0,imax+1,0,jmax+1,rho*ux);
  init_matrix(YMOM,0,imax+1,0,jmax+1,rho*uy);
  init_matrix(ENER,0,imax+1,0,jmax+1,0.5*rho*(ux*ux+uy*uy)+p/(k-1));
}

void init_sw(
  double h,
  double ux,
  double uy,
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM
)
{
  init_matrix(MASS,0,imax+1,0,jmax+1,h);
  init_matrix(XMOM,0,imax+1,0,jmax+1,h*ux);
  init_matrix(YMOM,0,imax+1,0,jmax+1,h*uy);
}


void init_flag_sw(
  char* geometry,
  int imax,
  int jmax,
  int **flag
)
{
   int xsize = imax + 2;
   int ysize = jmax + 2;
  int**	pic = read_pgm(geometry);

	for (int i=0; i<xsize; i++){
		for (int j=0; j<ysize; j++){

		    flag[i][j] = 0;

    		switch(pic[i][j]){
          //fluid
    			case 0: flag[i][j] = 1<<1; break;
          //no-slip
    			case 1: flag[i][j] = 1<<2; break;
          //free-slip
    			case 2: {
          flag[i][j]  = 1<<3;
          flag[i][j] |= FLUID; break;}
          //outflow
    			case 3: flag[i][j] = 1<<4; break;
          //inflow
    			case 4: flag[i][j] = 1<<0; break;
  		}
		}
	}

  for(int i = 0;i<xsize;++i){
    for(int j = 0;j<ysize;++j){
      if(flag[i][j] & NO_SLIP){
          if(i != 0       && (flag[i-1][j] & FLUID)) flag[i][j] |= B_W;
          if(i != xsize-1 && (flag[i+1][j] & FLUID)) flag[i][j] |= B_O;
          if(j != 0       && (flag[i][j-1] & FLUID)) flag[i][j] |= B_S;
          if(j != ysize-1 && (flag[i][j+1] & FLUID)) flag[i][j] |= B_N;
      }
    }
  }




	free_imatrix(pic, 0,xsize-1,0,ysize-1);

}

void init_flag(
  char* geometry,
  int imax,
  int jmax,
  int **flag
)
{
   int xsize = imax + 2;
   int ysize = jmax + 2;
  int**	pic = read_pgm(geometry);

	for (int i=0; i<xsize; i++){
		for (int j=0; j<ysize; j++){

		    flag[i][j] = 0;

    		switch(pic[i][j]){
          //fluid
    			case 0: flag[i][j] = 1<<1; break;
          //no-slip
    			case 1: flag[i][j] = 1<<2; break;
          //free-slip
    			case 2: flag[i][j] = 1<<3; break;
          //outflow
    			case 3: flag[i][j] = 1<<4; break;
          //inflow
    			case 4: flag[i][j] = 1<<0; break;
  		}
		}
	}

  for(int i = 0;i<xsize;++i){
    for(int j = 0;j<ysize;++j){
      if(flag[i][j] & NO_SLIP){
          if(i != 0       && (flag[i-1][j] & FLUID)) flag[i][j] |= B_W;
          if(i != xsize-1 && (flag[i+1][j] & FLUID)) flag[i][j] |= B_O;
          if(j != 0       && (flag[i][j-1] & FLUID)) flag[i][j] |= B_S;
          if(j != ysize-1 && (flag[i][j+1] & FLUID)) flag[i][j] |= B_N;
      }
    }
  }




	free_imatrix(pic, 0,xsize-1,0,ysize-1);

}

// void init_flag_sw(
//   char* geometry,
//   int imax,
//   int jmax,
//   int **flag
// )
// {
//    int xsize = imax + 2;
//    int ysize = jmax + 2;
//   int**	pic = read_pgm(geometry);
//
// 	for (int i=0; i<xsize; i++){
// 		for (int j=0; j<ysize; j++){
//
// 		    flag[i][j] = 0;
//
//     		switch(pic[i][j]){
//           //fluid
//     			case 0: flag[i][j] = 1<<1; break;
//           //no-slip
//     			case 1: flag[i][j] = 1<<2; break;
//           //free-slip
//     			case 2: {
//           flag[i][j]  = 1<<3;
//           flag[i][j] |= FLUID; break;}
//           //outflow
//     			case 3: flag[i][j] = 1<<4; break;
//           //inflow
//     			case 4: flag[i][j] = 1<<0; break;
//   		}
// 		}
// 	}
//
//   for(int i = 0;i<xsize;++i){
//     for(int j = 0;j<ysize;++j){
//       if(flag[i][j] & NO_SLIP){
//           if(i != 0       && (flag[i-1][j] & FLUID)) flag[i][j] |= B_W;
//           if(i != xsize-1 && (flag[i+1][j] & FLUID)) flag[i][j] |= B_O;
//           if(j != 0       && (flag[i][j-1] & FLUID)) flag[i][j] |= B_S;
//           if(j != ysize-1 && (flag[i][j+1] & FLUID)) flag[i][j] |= B_N;
//       }
//     }
//   }
//
//
//
//
// 	free_imatrix(pic, 0,xsize-1,0,ysize-1);
//
// }
