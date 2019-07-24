#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
#include <stdio.h>
#include "mpi.h"
#include"parallel.h"
#include <string.h>



int main(int argc, char* argv[])
{
    const char* filename= "cavity100.dat";

    //initialize variables
  	double t = 0; /*time start*/
  	int    n = 0; /*iteration and time step counter*/
    int    it;
  	double res; /*residual for SOR*/
    int time_count = 0;
  	/*arrays*/
  	double **U, **V, **P;
  	double **RS, **F, **G;
  	/*those to be read in from the input file*/
  	double Re, t_end, dt, dx, dy, alpha, omg, tau, eps, dt_value;
    double UI,VI,PI,GX,GY;
    int itermax,imax,jmax;
    double xlength, ylength;
    int  ir,il, jb, jt,omg_i, omg_j;
    int myrank,nproc;
    int rank_b,rank_l,rank_r,rank_t;
    int iproc,jproc;

  	/*MPI init*/
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &nproc ); /* asking for the number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank); /* asking for the local process id */
    MPI_Status status;
    /*read data from main process then broadcast them to all the other*/
    if(myrank == 0){
      printf("\nWelcome to Computational Fluid Dynamics Parallelization section!\n");
      printf("In this worksheet we will implement a distributed memory parallelization for the Navier-Stokes Equations\n");
      printf("Your visualization files are saved in the results folder\n");
      printf("Let us start!\n");
      read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax,
      &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
      if(iproc*jproc!=nproc){ iproc = nproc;jproc = 1;}
    }
    MPI_Bcast( &Re      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &UI      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &VI      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &PI      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &GX      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &GY      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &t_end   , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &xlength , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &ylength , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dt      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dx      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dy      , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &imax    , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &jmax    , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &alpha   , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &omg     , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &tau     , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &itermax , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &eps     , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dt_value, 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &jproc   , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &iproc   , 1,MPI_DOUBLE, 0, MPI_COMM_WORLD );

    /*initializationg of uvp*/
    init_parallel( iproc, jproc, imax, jmax, &myrank,&il,&ir,&jb,&jt,&rank_l,&rank_r,&rank_b,&rank_t,&omg_i,&omg_j,nproc);
    //printf("myrank %d il %d ir %d jb %d jt %d  ",myrank,il,ir,jb,jt);
    //printf("rank_l %d rank_r %d rank_b %d rank_t %d ",rank_l,rank_r,rank_b,rank_t);
    //printf("omg_i %d omg_j %d \n",omg_i,omg_j);
    /*allocate each small block*/

    P  = matrix(il-1, ir+1, jb-1, jt+1);
  	U  = matrix(il-2, ir+1, jb-1, jt+1);
  	V  = matrix(il-1, ir+1, jb-2, jt+1);
  	F  = matrix(il-2, ir+1, jb-1, jt+1);
  	G  = matrix(il-1, ir+1, jb-2, jt+1);
  	RS = matrix(il  , ir  , jb  , jt  );
   // init_uvp(UI, VI, PI,ir,il,jt,jb,U,V,P);
   // init_uvp(myrank,myrank,myrank,ir,il,jt,jb,U,V,P);

   /*initialize Buffers*/
   int bSize = max(omg_i+1,omg_j+1);
   double bufSend[bSize];
   double bufRecv[bSize];

   init_uvp(UI, VI, PI,ir,il,jt,jb,U,V,P);
   //printf("ID = %d has dt = %f\n", myrank, dt);
   //going through all time steps
   while(t < t_end){
        /*print running time*/
        if(myrank==0){
        time_count++;
        if ((int)(time_count*dt_value)%100==0){printf("The program is writing result for %f seconds \n",(time_count*dt_value/100));}
        }

        //adaptive time stepping
        double uMaxLocal = mmax(U,ir,il,jb,jt);
        double vMaxLocal = mmax(V,ir,il,jb,jt);
        double uMaxGlobal;
        double vMaxGlobal;
        MPI_Allreduce(&uMaxLocal,&uMaxGlobal,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&vMaxLocal,&vMaxGlobal,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        calculate_dt(Re,tau,&dt,dx,dy,uMaxGlobal,vMaxGlobal);

        //setting bound.values
        boundaryvalues(ir,il,jt,jb,rank_r,rank_l,rank_t,rank_b, U, V);

        //computing F, G and right hand side of pressue equation
        calculate_fg(Re, GX, GY, alpha, dt, dx, dy, ir,il,jt,jb,rank_r,rank_l,rank_t,rank_b, U, V, F, G);
        calculate_rs(dt, dx, dy,ir,il,jt,jb, F, G, RS);

        //iteration counter
        it = 0;
        do{
              //perform SOR iteration, at same time set bound.values for P and new residual value
              sor(omg, dx, dy,ir,il,jt,jb,rank_r,rank_l,rank_t,rank_b,omg_i,omg_j, P, RS, &res);
              pressure_comm(P, il,ir,jb,jt,rank_l,rank_r,rank_b,rank_t,bufSend,bufRecv, &status, myrank);
              it++;
              double sendRes = res;
              int sendIt  = it;
              MPI_Allreduce(&sendRes,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
              MPI_Allreduce(&sendIt,&it,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
         }while(it<itermax && res>eps);

      //calculate U and V of this time step
      calculate_uv(dt, dx, dy,ir,il,jt,jb, U, V, F, G, P);
      uv_comm(U,V,il,ir,jb,jt,rank_l,rank_r,rank_b,rank_t,bufSend,bufRecv,&status,myrank);
      t += dt;
      n++;
      //output of pics for animation
   }
  /*assign output file*/
  char sss[10];
  sprintf(sss, "%01d", myrank);
  char *output_file =  concat("results/Solution_",sss);
  output_uvp(U,V,P,n,dx,dy,il,ir,jb,jt,omg_i,omg_j,output_file);

  /*deallocate memory*/
  free_matrix(U , il-2, ir+1, jb-1, jt+1);
  free_matrix(V , il-1, ir+1, jb-2, jt+1);
  free_matrix(P , il-1, ir+1, jb-1, jt+1);
  free_matrix(RS, il  , ir  , jb  , jt  );
  free_matrix(F , il-2, ir+1, jb-1, jt+1);
  free_matrix(G , il-1, ir+1, jb-2, jt+1);

  MPI_Finalize();
  return 0;
}
