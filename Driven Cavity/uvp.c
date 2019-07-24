#include"helper.h"
#include"uvp.h"
#include<math.h>


void calculate_fg(
  double Re, // rexnolds Number
  double GX, // gravity x
  double GY, //gravity y
  double alpha, // constant
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
  ){

  double uij, vij, ulj, urj, uiu, uid, viu, vid, vrj, vlj, ulu, vrd;
	double Dx = 1/dx, Dy = 1/dy;
	int i, j;
	for (i=1; i<imax+1; i++){
		for  (j=1; j<jmax+1; j++){
			uij = U[i][j];
			ulj = U[i+1][j];
			urj = U[i-1][j];
			uiu = U[i][j+1];
			uid = U[i][j-1];
			ulu = U[i-1][j+1];

			vij = V[i][j];
			viu = V[i][j+1];
			vid = V[i][j-1];
			vrj = V[i+1][j];
			vlj = V[i-1][j];
			vrd = V[i+1][j-1];

			F[i][j] = uij + dt*(1/Re*((ulj-2*uij+urj)*Dx*Dx + (uiu-2*uij+uid)*Dy*Dy) -
				0.25*Dx*(pow((uij+ulj),2)-pow((urj+uij),2) + alpha*(fabs(uij+ulj)*(uij-ulj)-fabs(urj+uij)*(urj-uij))) -
				0.25*Dy*((vij+vrj)*(uij+uiu)-(vid+vrd)*(uid+uij) + alpha*(fabs(vij+vrj)*(uij-uiu)-fabs(vid+vrd)*(uid-uij))) +
				GX);

			G[i][j] = vij + dt*(1/Re*((vrj-2*vij+vlj)*Dx*Dx + (viu-2*vij+vid)*Dy*Dy) -
				0.25*Dy*(pow((vij+viu),2)-pow((vid+vij),2) + alpha*(fabs(vij+viu)*(vij-viu)-fabs(vid+vij)*(vid-vij))) -
				0.25*Dx*((uij+uiu)*(vij+vrj)-(urj+ulu)*(vlj+vij) + alpha*(fabs(uij+uiu)*(vij-vrj)-fabs(urj+ulu)*(vlj-vij))) +
				GY);

		}
		/* rewrite G(i,0) and G(i, jmax) with bound.cond. for G */
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
	for (j=1; j<jmax+1; j++){
		/* rewrite F(0,j) and F(imax, j) with bound.cond. for F */
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
} // end of funtion


  void calculate_rs(
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **F,
    double **G,
    double **RS
  ){
    double Dt = 1/dt;
    for(int i = 1; i<=imax; ++i){
      for( int j = 1; j<=jmax; ++j){
        RS[i][j] = Dt*((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy);
      } // end of j-loop
    } // end of i-loop
  } // end of function


  void calculate_dt(
    double Re,
    double tau,
    double *dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **U,
    double **V
  ){
  	double tmp, maxi1, maxi2;
  	/* we rewrite dt only if tau is nonnegative, else we do nothing */
  	if(tau>=0){
  		tmp = 0.5*Re*(dx*dx*dy*dy)/(dx*dx+dy*dy);
  		maxi1 = mmax(U, imax, jmax);
  		maxi2 = mmax(V, imax, jmax);
  		*dt = tau * fmin(tmp, fmin(dx/maxi1, dy/maxi2)); // (dx/mmax(U, imax, jmax), dy/mmax(V, imax, jmax));
   	}
  }


  void calculate_uv(
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    double **U,
    double **V,
    double **F,
    double **G,
    double **P
  ){
    for(int i = 1; i <= imax-1; ++i){
      for(int j = 1; j <= jmax; ++j){
        U[i][j] = F[i][j]-(dt/dx)*(P[i+1][j]-P[i][j]);

      }
    }
    for(int i = 1; i <= imax; ++i){
      for(int j = 1; j <= jmax-1; ++j){
          V[i][j]=G[i][j] - (dt/dy)*(P[i][j+1]-P[i][j]);
      }
    }
  } // end of function.
