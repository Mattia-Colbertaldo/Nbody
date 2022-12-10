#include "common.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <mpi.h>

int rank;

#define OK std::cout << "At main:" << __LINE__ << " from process " << rank << std::endl

// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream& fsave, std::vector<particle_pos>& parts, int num_parts, double size) {
    //int num_parts = parts.size();
    static bool first = true;

    if (first) {
        fsave << num_parts << " " << size << "\n";
        first = false;
    }

    for (int i = 0; i < num_parts; ++i) {
        fsave << parts[i].x << " " << parts[i].y << "\n";
    }

    fsave << std::endl;
}

// Particle Initialization
void init_particles(std::vector<particle_vel_acc>& parts_vel_acc_loc , std::vector<float>& masses, int num_parts, double size,int part_seed, int rank, std::vector<particle_pos>&parts_pos, int num_loc) {
    
    //int num_parts = parts.size();
    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());
    OK;
    int sx = (int)ceil(sqrt((double)num_parts));
    int sy = (num_parts + sx - 1) / sx;
    OK;

    std::vector<int> shuffle(num_parts);
    for (int i = 0; i < shuffle.size(); ++i) {
        shuffle[i] = i;
    }

    // initialize local vector of positions and velocities (parts_pos_vel_loc)
    // also fill the local part of vector of positions (parts_pos) and then allgather it
    // --> result : all have local values of positions and velocities in parts_pos_vel_loc and the positions of ALL particles in parts_pos
    for (int i = 0; i < num_loc; ++i) {
        // Make sure particles are not spatially sorted
        std::uniform_int_distribution<int> rand_int(0, num_parts - i - 1);
        int j = rand_int(gen);
        int k = shuffle[j];
        shuffle[j] = shuffle[num_parts - i - 1];

        // Distribute particles evenly to ensure proper spacing
        
        parts_pos[i+ num_loc].x = size * (1. + (k % sx)) / (1 + sx);
        parts_pos[i+ num_loc].y = size * (1. + (k / sx)) / (1 + sy);
        

        
        // Assign random velocities within a bound
        std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
        parts_vel_acc_loc[i].vx = rand_real(gen);
        parts_vel_acc_loc[i].vy = rand_real(gen);
        
        
    }

    //if(rank==0){
    //just process 0 init masses
    for(int i=0; i<num_parts; i++){
        std::uniform_real_distribution<float> rand_mass(0.001, 0.1);
        float m = rand_mass(gen);
        masses[i]=m;
    //    }
    }
    OK;
    MPI_Bcast(&masses , num_parts , MPI_FLOAT , 0 , MPI_COMM_WORLD); //FLAG BCAST MASSES
    MPI_Bcast(&parts_pos , num_parts , MPI_FLOAT , 0 , MPI_COMM_WORLD);

    MPI_Barrier( MPI_COMM_WORLD);
    OK;
    
    
    MPI_Allgather( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL ,  &parts_pos[0] , 2*num_loc , MPI_DOUBLE , MPI_COMM_WORLD);

}





// Command Line Option Processing
int find_arg_idx(int argc, char** argv, const char* option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::atoi(argv[iplace + 1]);
    }

    return default_value;
}

char* find_string_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}

// ==============
// Main Function
// ==============

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    //int rank;
    int mpi_size;
    
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int part_seed;
    double size;
    int num_parts;
    int num_th;
    char* savename;

    if(rank==0){
    
    // Open Output File
        savename = find_string_option(argc, argv, "-o", nullptr);
        if (savename != nullptr) std::cout << "Creating file " << savename << "..." << std::endl;
        //std::ofstream fsave(savename);
        if (savename != nullptr) std::cout << "File created." << std::endl;
       

        // Parse Args
        if (find_arg_idx(argc, argv, "-h") >= 0) {
            std::cout << "Options:" << std::endl;
            std::cout << "-h: see this help" << std::endl;
            std::cout << "-n <int>: set number of particles" << std::endl;
            std::cout << "-o <filename>: set the output file name" << std::endl;
            std::cout << "-s <int>: set particle initialization seed" << std::endl;
            std::cout << "-t <int>: set number of threads (working only in parallel mode) [default = 8]" << std::endl;
            return 0;
        }

    
        // Initialize Particles
        num_parts = find_int_arg(argc, argv, "-n", 1000);
        part_seed = find_int_arg(argc, argv, "-s", 0);
        size = sqrt(density * num_parts);
        num_th = find_int_arg(argc, argv, "-t", 8);
        
        // Assing random mass
        
        //std::cout << "mass: " <<  masses[i] << std::endl;

    
   }

    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `rank`
    // in this way we split the load in the most equilibrate way
    
   
    int num_loc = num_parts/mpi_size;
    std::vector<particle_vel_acc> parts_vel_acc_loc(num_loc);
    std::vector<particle_pos> parts_pos(num_parts);
    
    std::vector<float> masses(num_parts);
    
    std::cout << "Trying to init particles..." << std::endl;
    init_particles(parts_vel_acc_loc, masses, num_parts, size, part_seed, rank, parts_pos, num_loc); 

    
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();
    
  

    std::cout << "Trying to init simulation..." << std::endl;
    if(rank==0) init_simulation(parts_pos, masses, num_parts, size);
    MPI_Barrier( MPI_COMM_WORLD);
    std::cout << "Init simulation ended." << std::endl;

    auto init_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff_1 = init_time - start_time;
    double seconds_1 = diff_1.count();
    std::cout << "initialization Time = " << seconds_1 << " seconds\n";



    //for nel tempo: non parallelizzare
    for (int step = 0; step < nsteps; ++step) {
    
        
        simulate_one_step(parts_pos, parts_vel_acc_loc, masses, num_parts, num_loc, size);
        MPI_Barrier( MPI_COMM_WORLD);
        
        // Save state if necessary
        if(rank==0)
        {
            std::ofstream fsave(savename);
            if (fsave.good() && (step % savefreq) == 0) {
                //gather_for_save(parts, masses, num_parts, rank, size );
                save(fsave, parts_pos, num_parts, size);
            }
            
            if(step > 0){
                if (step%10 == 0){
                fflush(stdout);
                printf("[ %d% ]\r", (int)(step*100/nsteps));
                if(step = nsteps-1){
                    fsave.close();
                }
                }
            }
        }
    
        
    }
        
    

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();
    // Finalize
    //fsave.close();
    std::cout << "Simulation Time = " << seconds << " seconds for " << num_parts <<
     " particles and " << nsteps << " steps.\n";
    //delete[] parts;
    MPI_Finalize();
}
