#include"helper.h"
#include"boundary_val.h"
#include "mpi.h"

void boundaryvalues(
  int ir,
  int il,
  int jt,
  int jb,
  int rank_r,
  int rank_l,
  int rank_t,
  int rank_b,
  double **U,
  double **V)
  {
//U=[il-2,ir+1]×[jb-1,jt+1] V=[il-1,ir+1]×[jb-2,jt+1]
if(rank_l == MPI_PROC_NULL)
  for(int j=jb; j<=jt; ++j){
    //Equation 15
    U[il-2][j] = 0.0;
    V[il-1][j] = -V[il][j];
  }

if(rank_r == MPI_PROC_NULL)
    for(int j=jb; j<=jt; ++j){
      V[ir+1][j] = -V[ir][j];
      U[ir+1][j]    = 0.0;
    }
//
if(rank_b == MPI_PROC_NULL)
  for(int i = il; i<=ir; ++i){
    V[i][jb-2] = 0.0;
    U[i][jb-1] = -U[i][jb];
}

if(rank_t == MPI_PROC_NULL)
  for(int i = il; i<=ir; ++i){
    V[i][jt+1] = 0.0;
    U[i][jt+1] = 2.0 - U[i][jt]; // moving wall---here would be 2 times U_wall.
  }

}//end of function.
