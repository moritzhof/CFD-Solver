#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#include <iostream>


void init_parallel(int iproc,  int jproc,  int imax, int jmax,
                   int myrank, int il,     int ir, int jb, int jt,
                   int rank_l, int rank_r, int rank_b, int rank_t,
                   int omg_i,  int omg_j,  int num_proc){

      int myRank = myrank;

      
