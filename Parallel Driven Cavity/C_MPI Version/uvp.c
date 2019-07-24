#include"helper.h"
#include"uvp.h"
#include<math.h>
#include "mpi.h"


void calculate_fg(
  double Re, // rexnolds Number
  double GX, // gravity x
  double GY, //gravity y
  double alpha, // constant
  double dt,
  double dx,
  double dy,
  int ir,
  int il,
  int jt,
  int jb,
  int rank_r,
  int rank_l,
  int rank_t,
  int rank_b,
  double **U,
  double **V,
  double **F,
  double **G
  ){
//F=[il-2,ir+1]×[jb-1,jt+1]  G=[il-1,ir+1]×[jb-2,jt+1]
  double uij, vij, ulj, urj, uiu, uid, viu, vid, vrj, vlj, ulu, vrd;
	double Dx = 1/dx, Dy = 1/dy;
	int i, j;
	for (i=il-1; i<=ir; i++){
		for  (j=jb; j<=jt; j++){
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
			vrd = V[i+1][j-1];

			F[i][j] = uij + dt*(1/Re*((ulj-2*uij+urj)*Dx*Dx + (uiu-2*uij+uid)*Dy*Dy) -
				0.25*Dx*(pow((uij+ulj),2)-pow((urj+uij),2) + alpha*(fabs(uij+ulj)*(uij-ulj)-fabs(urj+uij)*(urj-uij))) -
				0.25*Dy*((vij+vrj)*(uij+uiu)-(vid+vrd)*(uid+uij) + alpha*(fabs(vij+vrj)*(uij-uiu)-fabs(vid+vrd)*(uid-uij))) +
				GX);
      }
		}

    for (i=il; i<=ir; i++){
      for  (j=jb-1; j<=jt; j++){
        uij = U[i][j];
        ulj = U[i+1][j];
        urj = U[i-1][j];
        uiu = U[i][j+1];
        ulu = U[i-1][j+1];

        vij = V[i][j];
        viu = V[i][j+1];
        vid = V[i][j-1];
        vrj = V[i+1][j];
        vlj = V[i-1][j];
        vrd = V[i+1][j-1];


        G[i][j] = vij + dt*(1/Re*((vrj-2*vij+vlj)*Dx*Dx + (viu-2*vij+vid)*Dy*Dy) -
          0.25*Dy*(pow((vij+viu),2)-pow((vid+vij),2) + alpha*(fabs(vij+viu)*(vij-viu)-fabs(vid+vij)*(vid-vij))) -
          0.25*Dx*((uij+uiu)*(vij+vrj)-(urj+ulu)*(vlj+vij) + alpha*(fabs(uij+uiu)*(vij-vrj)-fabs(urj+ulu)*(vlj-vij))) +
          GY);

      }
		/* rewrite G(i,0) and G(i, jmax) with bound.cond. for G */

	}

if(rank_b == MPI_PROC_NULL)
for (i=il; i<=ir; i++){
  G[i][jb-1] = V[i][jb-1];
}

if(rank_t == MPI_PROC_NULL)
for (i=il; i<=ir; i++){
  G[i][jt] = V[i][jt];
}
if(rank_l == MPI_PROC_NULL)
	for (j=jb; j<=jt; j++){
		/* rewrite F(0,j) and F(imax, j) with bound.cond. for F */
		F[ir-1][j] = U[ir-1][j];
  }
if(rank_r == MPI_PROC_NULL)
	for (j=jb; j<=jt; j++){
		F[ir][j] = U[ir][j];
	}
} // end of funtion


  void calculate_rs(
    double dt,
    double dx,
    double dy,
    int ir,
    int il,
    int jt,
    int jb,
    double **F,
    double **G,
    double **RS
  ){
    //RS=[il,ir]×[jb,jt]
    double Dt = 1/dt;
    for(int i = il; i<=ir; ++i){
      for( int j = jb; j<=jt; ++j){
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
    double uMax,
    double vMax
  ){
    double tmp;
  	/* we rewrite dt only if tau is nonnegative, else we do nothing */
  	if(tau>=0){
  		tmp = 0.5*Re*(dx*dx*dy*dy)/(dx*dx+dy*dy);
  		*dt = tau * fmin(tmp, fmin(dx/uMax, dy/vMax)); // (dx/mmax(U, imax, jmax), dy/mmax(V, imax, jmax));
   	}
  }

//index might be wrong
  void calculate_uv(
    double dt,
    double dx,
    double dy,
    int ir,
    int il,
    int jt,
    int jb,
    double **U,
    double **V,
    double **F,
    double **G,
    double **P
  ){
    //U=[il-2,ir+1]×[jb-1,jt+1] V=[il-1,ir+1]×[jb-2,jt+1]
    for(int i = il-1; i <= ir; ++i){
      for(int j = jb; j <= jt; ++j){
        U[i][j] = F[i][j]-(dt/dx)*(P[i+1][j]-P[i][j]);

      }
    }
    for(int i = il; i <= ir; ++i){
      for(int j = jb-1; j <= jt; ++j){
          V[i][j]=G[i][j] - (dt/dy)*(P[i][j+1]-P[i][j]);
      }
    }
  } // end of function.
