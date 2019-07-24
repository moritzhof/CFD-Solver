#include <stdio.h>
#include <mpi.h>

void init_parallel(int iproc,int jproc,int imax,int jmax,
int *myrank,int *il,int *ir,int *jb,int *jt,
int *rank_l,int *rank_r,int *rank_b,int *rank_t,
int *omg_i,int *omg_j,int num_proc){
  printf("%d",*myrank);
}


int main(int argc, char* argv[]){
  MPI_Init( &argc, &argv ); /* execute n processes */
  int iproc  = 1;
  int jproc = 2;
  int *il,*ir,*jt,*jb,*omg_i,*omg_j;
  int *rank_b,*rank_l,*rank_t,*rank_r;
//  int myrank, nproc;
int* nproc;
int* myrank;



  int imax = 103;
  int jmax = 50;

  MPI_Comm_rank( MPI_COMM_WORLD, myrank); /* asking for the local process id */

  MPI_Comm_size( MPI_COMM_WORLD, nproc); /* asking for the number of processes */


void init_parallel(iproc,jproc,imax,jmax, temp,il,ir,jb,jt,rank_l,rank_r,rank_b,rank_t,omg_i,omg_j,nproc);

  MPI_Finalize();
}
