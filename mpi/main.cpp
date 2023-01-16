#include "common.h"

#include "Find_Arg.hpp"
#include "PhysicalForce.hpp"
#include "Output.hpp"
#include "Particle.hpp"
#include "Simulation.hpp"

#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <mpi.h>


using namespace common_h;

int rank, mpi_size;
bool first = true;



#define OK std::cout << "At main:" << __LINE__ << " from process " << rank << std::endl




// ==============
// Main Function
// ==============

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    //int rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //std::cout << "mpi_size: " << mpi_size << std::endl;

    int part_seed;
    double size;
    int num_parts;
    int num_th;

    Find_Arg finder= Find_Arg(argc, argv);
    std::string savename = finder.find_string_option("-o", "out.txt");
    if(rank==0 && savename != "") std::cout << "Creating file " << savename << "..." << std::endl;
    std::ofstream fsave(savename);
    if (rank == 0 && savename != "") std::cout << "File created." << std::endl;
    if(rank != 0) fsave.close();

    std::vector<int> sizes(mpi_size);
    std::vector<int> displs(mpi_size+1);
    
    std::string forcename;
    
    if(rank==0){       
    
    // Open Output File
        finder.find_int_arg("-h", 0);

        // Initialize Particles
        num_parts = finder.find_int_arg("-n", 1000);
    
        const int part_seed = finder.find_int_arg("-s", 0);
        size = std::sqrt(density * num_parts);
        const int num_th = finder.find_int_arg("-t", 8);

        displs[0]=0;
        for(int i=0; i<mpi_size; i++){
            sizes[i]= num_parts/mpi_size + (i < num_parts%mpi_size);
            //std::cout << "Size of process " << i << " = " << sizes[i] << std::endl;
            displs[i+1]= displs[i]+sizes[i];
            //std::cout << "Displ of process " << i << " = " << displs[i] << std::endl;
            
        }
   }

   

    //Find force
    forcename = finder.find_string_option("-f", "repulsive");

    std::unique_ptr<AbstractForce> force= finder.find_force(forcename);
    
    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `rank`
    // in this way we split the load in the most equilibrate way
    int num_loc, displ_loc;
    MPI_Bcast( &num_parts , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
   
    
    MPI_Bcast( &size , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &sizes[0] , mpi_size , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &displs[0] , mpi_size+1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &part_seed , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    num_loc=sizes[rank];
    displ_loc=displs[rank];
   
    Simulation simulation = Simulation(num_parts, num_loc);
    
    if (!rank) std::cout << "Initialization: ";
    MPI_Datatype mpi_parts_pos_type=simulation.init_particles(num_parts, size, part_seed, num_loc, displs, sizes, rank);

    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

   if(!rank){
        auto init_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff_1 = init_time - start_time;
        double seconds_1 = diff_1.count();
        std::cout << seconds_1 << " seconds\n";
    }

    Output output= Output(savename);
  
    if(rank == 0) output.save(simulation.parts_pos , size, nsteps);

    //for nel tempo: non parallelizzare
    for (int step = 0; step < nsteps; ++step) {
        simulation.simulate_one_step(num_parts, num_loc, displ_loc, size, force);
        
        MPI_Barrier( MPI_COMM_WORLD);

        // Allgather delle posizioni, in questo modo aggiorno la posizione di tutte le particelle per tutti i processori.
        // Non serve comunicare velocità e accelerazione visto che sono necessarie solo localmente.
        MPI_Allgatherv( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL , &simulation.parts_pos[0] , &sizes[0] , &displs[0] , mpi_parts_pos_type ,
                        MPI_COMM_WORLD);
        
        // Save state if necessary
        if(rank==0)
        {
            output.save_output(savefreq, simulation.parts_pos , step, nsteps, size);
        }
        
    }

    MPI_Barrier( MPI_COMM_WORLD);

    if(!rank){
        auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    // Finalize

    std::cout << "Closing file..." << std::endl;
    fsave.close();

    
    std::cout << "Simulation Time = " << seconds << " seconds for " << num_parts <<
     " particles and " << nsteps << " steps.\n";
    //delete[] parts;
    }
    MPI_Finalize();
    

    
}
