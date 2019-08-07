// #include <boost/mpi/environment.hpp>
// #include <boost/mpi/communicator.hpp>
// #include <boost/mpi.hpp>

#include <iostream>
#include <mpi.h>
#include <vector>

    typedef std::vector<double> Vector;
    typedef std::vector<Vector> Matrix;

    typedef std::vector<int> intVector;
    typedef std::vector<intVector> intMatrix;


int main(int argc, char* argv[]){

    //--------------------------------------------
    //            size of the domain
    //--------------------------------------------
    int xlength = 1;
    int ylength	= 1;

    // #--------------------------------------------
    // #            number of cells
    // #--------------------------------------------
    int imax = 50;
    int jmax = 50;

    // #--------------------------------------------
    // #               time steps
    // #--------------------------------------------
    double dt    = 0.05;
    double t_end = 25.0;
    double tau	 = 0.5;

    // #--------------------------------------------
    // #               output
    // #--------------------------------------------
    double dt_value = 0.5;

    // #--------------------------------------------
    // #               pressure
    // #--------------------------------------------
    int itermax  = 100;
    double eps	 = 0.001;
    double omg    = 1.7;
    double alpha = 0.5;

    // #--------------------------------------------
    // #               reynoldsnumber
    // #--------------------------------------------
    double Re = 100.0;

    // #--------------------------------------------
    // #               gravitation
    // #--------------------------------------------
    double GX = 0;
    double GY =	0;

    // #--------------------------------------------
    // #         initialization pressure
    // #--------------------------------------------
    double PI = 0;

    // #--------------------------------------------
    // #       initialization velocity
    // #--------------------------------------------
    double UI = 0;
    double VI =	0;

    // #--------------------------------------------
    // #         iproc & jproc
    // #--------------------------------------------
    int iproc = 2;
    int jproc = 1;
    int rank_b,rank_l,rank_r,rank_t;
    int dy, dx;

    // boost::mpi::environment env(argc, argv);
    // boost::mpi::communicator world;
    //
    // int myrank = world.rank();
    // int nproc  = world.size();

    MPI::Init(argc, argv);
    int myrank = MPI::COMM_WORLD.Get_rank();
    int nproc  = MPI::COMM_WORLD.Get_size();

    if(myrank == 0){
        std::cout<< " Welcome to Computational Fluid Dynamics Parallelization Solver "<<std::endl;
        std::cout<< " The Lid Driven Cavity is simulated" << std::endl;
        std::cout<< " Total Processes: " << nproc << std::endl;
    }

    MPI::COMM_WORLD.Bcast( &Re      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &UI      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &VI      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &PI      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &GX      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &GY      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &t_end   , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &xlength , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &ylength , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &dt      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &dx      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &dy      , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &imax    , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &jmax    , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &alpha   , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &omg     , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &tau     , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &itermax , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &eps     , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &dt_value, 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &jproc   , 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast( &iproc   , 1, MPI::DOUBLE, 0);


    MPI::Finalize();
    return 0;
}



    // boost MPI broadcast:
    // broadcast(world, Re, 0);
    // broadcast(world, UI, 0);
    // broadcast(world, VI, 0);
    // broadcast(world, PI, 0);
    // broadcast(world, GX, 0);
    // broadcast(world, GY, 0);
    // broadcast(world, t_end, 0);
    // broadcast(world, xlength, 0);
    // broadcast(world, ylength, 0);
    // broadcast(world, dt, 0);
    // broadcast(world, dx, 0);
    // broadcast(world, dy, 0);
    // broadcast(world, imax, 0);
    // broadcast(world, jmax, 0);
    // broadcast(world, alpha, 0);
    // broadcast(world, omg, 0);
    // broadcast(world, tau, 0);
    // broadcast(world, itermax, 0);
    // broadcast(world, eps, 0);
    // broadcast(world, dt_value, 0);
    // broadcast(world, iproc, 0);
    // broadcast(world, jproc, 0);
