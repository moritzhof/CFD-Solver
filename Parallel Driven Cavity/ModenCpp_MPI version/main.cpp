#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#include <iostream>


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

    boost::mpi::environment env;
    boost::mpi::communicator world;

    int myrank = world.rank();
    int nproc  = world.size();

    if(myrank == 0){
        std::cout<< " Welcome to Computational Fluid Dynamics Parallelization Solver "<<std::endl;
        std::cout<< " The Lid Driven Cavity is simulated" << std::endl;
        std::cout<< " Total Processes: " << nproc << std::endl;
    }

    boost::mpi::broadcast(world, Re, 0);


    MPI::Finalize();
    return 0;
}
