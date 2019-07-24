#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}



void init_parallel(int iproc,int jproc,int imax,int jmax,
int *myrank,int *il,int *ir,int *jb,int *jt,
int *rank_l,int *rank_r,int *rank_b,int *rank_t,
int *omg_i,int *omg_j,int num_proc){

  int myRank = *myrank;



  //Check Left Boundary
  if(myRank%iproc == 0){ *rank_l  = MPI_PROC_NULL;}
  else {*rank_l = myRank - 1;}

  //Check Right Boundary
  if(myRank%iproc == iproc - 1){ *rank_r = MPI_PROC_NULL; }
  else {*rank_r = myRank + 1;}


  //Check Top Boundary
  if((int)*myrank/iproc == 0){*rank_b = MPI_PROC_NULL;}
  else {*rank_b = *myrank - iproc;}
  //Check Bottom
  if((int)*myrank/iproc == jproc - 1) {*rank_t = MPI_PROC_NULL;}
  else {*rank_t = *myrank + iproc;}


  //
  if(*myrank == 0){
    if ( imax % iproc){*omg_i = imax / iproc + 1;}
    else { *omg_i = imax / iproc;}

    if ( jmax % jproc){*omg_j = jmax / jproc + 1;}
    else { *omg_j = jmax / jproc;}

    *il = 1;
    *ir = *omg_i;
    *jb = 1;
    *jt = *omg_j;
    for(int processID = 1; processID < num_proc; ++processID){

      int tempI;
      int tempJ;
      int IR;
      int IL;
      int JT;
      int JB;

      if ( processID % iproc < imax % iproc)
      { tempI = imax / iproc + 1;
        IL = tempI * (processID % iproc)+1;
        IR = IL + tempI-1;
      }
      else
      { tempI = imax / iproc;
        IL = imax % iproc * (tempI + 1) + (processID % iproc - imax % iproc) * tempI +1;
        IR = IL + tempI-1;
      }

      if ( processID / iproc < jmax % jproc)
      {tempJ = jmax / jproc + 1;
        JB = tempJ * (processID / iproc) + 1;
        JT = JB + tempJ-1;
      }
      else
      { tempJ = jmax / jproc;
        JB = jmax % jproc * (tempJ + 1) + (processID / iproc - jmax % jproc) * tempJ +1;
        JT = JB + tempJ-1;
      }

      MPI_Send(&tempI ,1,MPI_INT,processID,0,MPI_COMM_WORLD);
      MPI_Send(&tempJ,1,MPI_INT,processID,0,MPI_COMM_WORLD);
      MPI_Send(&IL ,1,MPI_INT,processID,0,MPI_COMM_WORLD);
      MPI_Send(&IR,1,MPI_INT,processID,0,MPI_COMM_WORLD);
      MPI_Send(&JT ,1,MPI_INT,processID,0,MPI_COMM_WORLD);
      MPI_Send(&JB,1,MPI_INT,processID,0,MPI_COMM_WORLD);

    }
  }
  else{
      MPI_Status status;
      MPI_Recv(omg_i,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
      MPI_Recv(omg_j,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
      MPI_Recv(il,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
      MPI_Recv(ir,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
      MPI_Recv(jt,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
      MPI_Recv(jb,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);

    }
}

void pressure_comm(
    double **P,
    int il,
    int ir,
    int jb,
    int jt,
    int rank_l,
    int rank_r,
    int rank_b,
    int rank_t,
    double *bufSend,
    double *bufRecv,
    MPI_Status *status,
    int chunk){
   int i;
    int j;

     //order might lead to problem
     int countj =jt-jb+1;
     int counti =ir-il+1;

      // send to the left | receive from the right,
      for(j = jb; j <= jt; j++) {bufSend[j-jb] = P[il][j];}
      MPI_Sendrecv(bufSend,countj,MPI_DOUBLE,rank_l,0,bufRecv,countj,MPI_DOUBLE,rank_r,0,MPI_COMM_WORLD,status);
      if(rank_r != MPI_PROC_NULL) for(j = jb; j <= jt; j++) { P[ir+1][j] = bufRecv[j-jb];}
      // send to the right | receive from the left,
      for(j = jb; j <= jt; j++) {bufSend[j-jb] = P[ir][j];}
      MPI_Sendrecv(bufSend,countj,MPI_DOUBLE,rank_r,0,bufRecv,countj,MPI_DOUBLE,rank_l,0,MPI_COMM_WORLD,status);
      if(rank_l != MPI_PROC_NULL) for(j = jb; j <= jt; j++) { P[il-1][j] = bufRecv[j-jb];}
      // send to the top | receive from the bottom,
      for(i = il; i <= ir; i++) {bufSend[i-il] = P[i][jt];}
      MPI_Sendrecv(bufSend,counti,MPI_DOUBLE,rank_t,0,bufRecv,counti,MPI_DOUBLE,rank_b,0,MPI_COMM_WORLD,status);
      if(rank_b != MPI_PROC_NULL) for(i = il; i <= ir; i++) {P[i][jb-1] = bufRecv[i -il];}
     // send to the bottom | receive from the top.
      for(i = il; i <= ir; i++) {bufSend[i-il] = P[i][jb];}
      MPI_Sendrecv(bufSend,counti,MPI_DOUBLE,rank_b,0,bufRecv,counti,MPI_DOUBLE,rank_t,0,MPI_COMM_WORLD,status);
      if(rank_t != MPI_PROC_NULL) for(i = il; i <= ir; i++) {P[i][jt+1] = bufRecv[i -il];}
}

void uv_comm(
  double**U,
  double**V,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double *bufSend,
  double *bufRecv,
  MPI_Status *status,
  int chunk){

    int j;
    int i;
    int countj_v = jt-jb+2;
    int countj_u = jt-jb+1;
    int counti_v = ir-il+1;
    int counti_u = ir-il+2;
    //U=[il-2,ir+1]×[jb-1,jt+1]; V=[il-1,ir+1]×[jb-2,jt+1];

    //send to the left | receive from the right,
    for(j = jb; j <= jt; j++) {bufSend[j-jb] = U[il-1][j];}
    MPI_Sendrecv(bufSend,countj_u,MPI_DOUBLE,rank_l,0,bufRecv,countj_u,MPI_DOUBLE,rank_r,0,MPI_COMM_WORLD,status);
    if(rank_r != MPI_PROC_NULL) for(j = jb; j <= jt; j++) { U[ir+1][j] = bufRecv[j-jb];}

    for(j = jb-2; j <= jt; j++) {bufSend[j-jb+1] = V[il][j];}
    MPI_Sendrecv(bufSend,countj_v,MPI_DOUBLE,rank_l,0,bufRecv,countj_v,MPI_DOUBLE,rank_r,0,MPI_COMM_WORLD,status);
    if(rank_r != MPI_PROC_NULL) for(j = jb-1; j <= jt; j++) { V[ir+1][j] = bufRecv[j-jb+1];}

    //send to the right | receive from the left,
    for(j = jb; j <= jt; j++) {bufSend[j-jb] = U[ir][j];}
    MPI_Sendrecv(bufSend,countj_u,MPI_DOUBLE,rank_r,0,bufRecv,countj_u,MPI_DOUBLE,rank_l,0,MPI_COMM_WORLD,status);
    if(rank_l != MPI_PROC_NULL) for(j = jb; j <= jt; j++) {U[il-2][j] = bufRecv[j-jb];}

    for(j = jb-1; j <= jt; j++) {bufSend[j-jb+1] = V[ir][j];}
    MPI_Sendrecv(bufSend,countj_v,MPI_DOUBLE,rank_r,0,bufRecv,countj_v,MPI_DOUBLE,rank_l,0,MPI_COMM_WORLD,status);
    if(rank_l != MPI_PROC_NULL) for(j = jb-1; j <= jt; j++) {V[il-1][j] = bufRecv[j-jb+1];}

    //send to the top | receive from the bottom
    for(i = il-1; i <= ir; i++) {bufSend[i-il+1] = U[i][jt];}
    MPI_Sendrecv(bufSend,counti_u,MPI_DOUBLE,rank_t,0,bufRecv,counti_u,MPI_DOUBLE,rank_b,0,MPI_COMM_WORLD,status);
    if(rank_b != MPI_PROC_NULL) for(i = il-1; i <= ir; i++) {U[i][jb-1] = bufRecv[i-il+1];}

    for(int i = il; i <= ir; i++) {bufSend[i-il] = V[i][jt];}
    MPI_Sendrecv(bufSend,counti_v,MPI_DOUBLE,rank_t,0,bufRecv,counti_v,MPI_DOUBLE,rank_b,0,MPI_COMM_WORLD,status);
    if(rank_b != MPI_PROC_NULL) for(i = il; i <= ir; i++) {V[i][jb-2] = bufRecv[i-il];}


   //send to the bottom | receive from the top.
    for(i = il-1; i <= ir; i++) {bufSend[i-il+1] = U[i][jb];}
    MPI_Sendrecv(bufSend,counti_u,MPI_DOUBLE,rank_b,0,bufRecv,counti_u,MPI_DOUBLE,rank_t,0,MPI_COMM_WORLD,status);
    if(rank_t != MPI_PROC_NULL) for(i = il-1; i <= ir;++i) {U[i][jt+1] = bufRecv[i-il+1];}

    for(i = il; i <= ir; i++) {bufSend[i-il] = V[i][jb-1];}
    MPI_Sendrecv(bufSend,counti_v,MPI_DOUBLE,rank_b,0,bufRecv,counti_v,MPI_DOUBLE,rank_t,0,MPI_COMM_WORLD,status);
    if(rank_t != MPI_PROC_NULL) for(i = il; i <= ir; i++) {V[i][jt+1] = bufRecv[i-il];}
}




  /*int MPI_Sendrecv(
  void *sendbuf,
  int sendcount,
  MPI_Datatype sendtype,
  int dest,
  int sendtag,
  void *recvbuf,
  int recvcount,
  MPI_Datatype recvtype,
  int source,
  int recvtag,
  MPI_Comm comm,
  MPI_Status *status
);*/
