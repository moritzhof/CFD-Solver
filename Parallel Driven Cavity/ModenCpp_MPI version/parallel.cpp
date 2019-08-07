
#include<iostream>
#include<mpi.h>

void init_parallel(int iproc,  int jproc,  int imax, int jmax,
                   int* myrank, int il,     int ir, int jb, int jt,
                   int rank_l, int rank_r, int rank_b, int rank_t,
                   int omg_i,  int omg_j,  int num_proc){

      int myRank = *myrank;
      //Check Left Boundary
      myRank%iproc == 0 ? rank_l = MPI::PROC_NULL : rank_l = myRank-1;
      //Check Right Boundary
      myRank%jproc == iproc-1 ? rank_r = MPI::PROC_NULL : rank_r = myRank+1;
      //Check Top Boundary
      (int)*myrank/iproc == 0 ? rank_b = MPI::PROC_NULL : rank_b = *myrank-iproc;
      //Check Bottem Boundary
      (int)*myrank/iproc == jproc-1 ? rank_t = MPI::PROC_NULL : rank_t = *myrank + iproc;
}

int main(int argc, char const *argv[]) {

  return 0;
}
