#include "sor.h"
#include <math.h>
#include "mpi.h"

void sor(
  double omg,
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
  int omg_i,
  int omg_j,
  double **P,
  double **RS,
  double *res
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration [il-1,ir+1]×[jb-1,jt+1]*/
  for(i = il; i <= ir; i++) {
    for(j = jb; j<=jt; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = il; i <= ir; i++) {
    for(j = jb; j <= jt; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
//  rloc = rloc/(imax*jmax);
rloc = rloc/(omg_i*omg_j);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;

  /* set boundary values [il-1,ir+1]×[jb-1,jt+1]*/
if(rank_l == MPI_PROC_NULL){
  for(j = jb; j <= jt; j++) {
    P[il-1][j] = P[il][j];
  }
}
if(rank_r == MPI_PROC_NULL){
  for(j = jb; j <= jt; j++){
    P[ir+1][j] = P[ir][j];
  }
}
if(rank_b == MPI_PROC_NULL){
  for(i = il; i <= ir; i++) {
    P[i][jb-1] = P[i][jb];
  }
}

if(rank_t == MPI_PROC_NULL){
    for(i = il; i <= ir; i++) {
      P[i][jt+1] = P[i][jt];
    }
  }



}
