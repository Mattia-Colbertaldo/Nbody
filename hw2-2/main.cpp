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

int rank, mpi_size;
bool first = true;
MPI_Datatype mpi_part_vel_acc_type;
MPI_Datatype mpi_part_pos_type;

#define OK std::cout << "At main:" << __LINE__ << " from process " << rank << std::endl

// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream& fsave, std::vector<particle_pos>& parts, int num_parts, double size) {
    //int num_parts = parts.size();

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
void init_particles(std::vector<particle_vel_acc>& parts_vel_acc_loc,
                    std::vector<float>& masses,
                    int num_parts, double size,int part_seed, int rank,
                    std::vector<particle_pos>&parts_pos, int num_loc, 
                    int loc_displ, std::vector<int>& sizes) {
    
    //int num_parts = parts.size();
    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());
    
    int sx = (int)ceil(sqrt((double)num_parts));
    int sy = (num_parts + sx - 1) / sx;


    ////// MPI STRUCT /////////////

    const int nitems=4;
    int          blocklengths[4] = {1,1,1,1};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    //MPI_Datatype mpi_part_vel_acc_type;
    MPI_Aint     offsets[4];

    offsets[0] = offsetof(particle_vel_acc, vx);
    offsets[1] = offsetof(particle_vel_acc, vy);
    offsets[2] = offsetof(particle_vel_acc, ax);
    offsets[3] = offsetof(particle_vel_acc, ay);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_part_vel_acc_type);
    MPI_Type_commit(&mpi_part_vel_acc_type);
    
    ///////////////////////////////

    ////// MPI STRUCT /////////////
    
    const int nitems1=2;
    int          blocklengths1[2] = {1,1};
    MPI_Datatype types1[2] = {MPI_DOUBLE, MPI_DOUBLE};
    //MPI_Datatype mpi_part_pos_type;
    MPI_Aint     offsets1[2];

    offsets1[0] = offsetof(particle_pos, x);
    offsets1[1] = offsetof(particle_pos, y);

    MPI_Type_create_struct(nitems1, blocklengths1, offsets1, types1, &mpi_part_pos_type);
    MPI_Type_commit(&mpi_part_pos_type);
    
    ///////////////////////////////


    std::vector<int> shuffle(num_parts);
    for (int i = 0; i < shuffle.size(); ++i) {
        shuffle[i] = i;
    }

    if(rank==0){

        std::vector<particle_vel_acc> parts_vel_acc_temp(num_parts);

        // initialize local vector of positions and velocities (parts_pos_vel_loc)
        // also fill the local part of vector of positions (parts_pos) and then allgather it
        // --> result : all have local values of positions and velocities in parts_pos_vel_loc and the positions of ALL particles in parts_pos
        for (int i = 0; i < num_parts; ++i) {
            // Make sure particles are not spatially sorted
            std::uniform_int_distribution<int> rand_int(0, num_parts - i - 1);
            int j = rand_int(gen);
            int k = shuffle[j];
            shuffle[j] = shuffle[num_parts - i - 1];
            // Distribute particles evenly to ensure proper spacing
            
            parts_pos[i].x = size * (1. + (k % sx)) / (1 + sx);
            parts_pos[i].y = size * (1. + (k / sx)) / (1 + sy);
            
            // Assign random velocities within a bound
            std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
            parts_vel_acc_temp[i].vx = rand_real(gen);
            parts_vel_acc_temp[i].vy = rand_real(gen);

            std::uniform_real_distribution<float> rand_mass(0.001, 0.1);
            float m = rand_mass(gen);
            masses[i]=m;

        }

        //Saving my data (i'm rank 0)
        for(int i=0; i<num_loc; i++){
            parts_vel_acc_loc[i].vx = parts_vel_acc_temp[i].vx;
            parts_vel_acc_loc[i].vy = parts_vel_acc_temp[i].vy;
        }
        //sending to others
        int d=sizes[0];
        for(int r=1; r<mpi_size; r++){
            MPI_Send( &parts_vel_acc_temp[d], sizes[r] , mpi_part_vel_acc_type , r , r, MPI_COMM_WORLD);
            d+=sizes[r];
        }
         
    }
    //else (i'm not rank 0)
    else MPI_Recv( &parts_vel_acc_loc[0] , num_loc , mpi_part_vel_acc_type , 0 , rank , MPI_COMM_WORLD , MPI_STATUS_IGNORE);


    
    
    //MPI_Barrier( MPI_COMM_WORLD );


    MPI_Bcast( &masses[0] , num_parts , MPI_FLOAT , 0 , MPI_COMM_WORLD);
    MPI_Bcast(&parts_pos[0] , num_parts, mpi_part_pos_type , 0 , MPI_COMM_WORLD);



    //MPI_Barrier( MPI_COMM_WORLD );
    
    
    //MPI_Allgather( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL ,  &parts_pos[0] , 2*num_loc, MPI_DOUBLE , MPI_COMM_WORLD);
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
    
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //std::cout << "mpi_size: " << mpi_size << std::endl;

    int part_seed;
    double size;
    int num_parts;
    int num_th;
    char* savename;

    savename = find_string_option(argc, argv, "-o", nullptr);
    if(rank==0 && savename != nullptr) std::cout << "Creating file " << savename << "..." << std::endl;
    std::ofstream fsave(savename);
    if (rank == 0 && savename != nullptr) std::cout << "File created." << std::endl;
    if(rank != 0) fsave.close();

    std::vector<int> sizes(mpi_size);
    std::vector<int> displs(mpi_size+1);
    
    

    if(rank==0){
        

    

    //savename = find_string_option(argc, argv, "-o", nullptr);
    //if (savename != nullptr) std::cout << "Creating file " << savename << "..." << std::endl;
    //std::ofstream fsave(savename);
    //if (savename != nullptr) std::cout << "File created." << std::endl;
       
    
    // Open Output File
        

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

        displs[0]=0;
        for(int i=0; i<mpi_size; i++){
            sizes[i]= num_parts/mpi_size + (i < num_parts%mpi_size);
            std::cout << "Size of process " << i << " = " << sizes[i] << std::endl;
            displs[i+1]= displs[i]+sizes[i];
            std::cout << "Displ of process " << i << " = " << displs[i] << std::endl;
            MPI_Send( &sizes[i] , 1 , MPI_INT , i , i+mpi_size , MPI_COMM_WORLD);
            MPI_Send( &displs[i] , 1 , MPI_INT , i , i , MPI_COMM_WORLD);
        }
        
        
    
   }

    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `rank`
    // in this way we split the load in the most equilibrate way
    int num_loc, displ_loc;
    MPI_Bcast( &num_parts , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &size , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &sizes[0] , mpi_size , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &displs[0] , mpi_size+1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &part_seed , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Recv( &num_loc , 1 , MPI_INT , 0 , rank+mpi_size , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
    MPI_Recv( &displ_loc , 1 , MPI_INT , 0 , rank , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
    
    //std::cout << "num_loc: " << num_loc << std::endl;
    //std::cout << "size: " << size << std::endl;
    std::vector<particle_vel_acc> parts_vel_acc_loc(num_loc);
    std::vector<particle_pos> parts_pos(num_parts);
    
    std::vector<float> masses(num_parts);
    
    if (!rank) std::cout << "Trying to init particles..." << std::endl;
    init_particles(parts_vel_acc_loc, masses, num_parts, size, part_seed, rank, parts_pos, num_loc, displ_loc, sizes); 
    
    
    // Algorithm
    //std::cout << rank << " -> Starting..." << std::endl;
    auto start_time = std::chrono::steady_clock::now();
    
    /*

    if(!rank) std::cout << "Trying to init simulation..." << std::endl;
    //init_simulation(parts_pos, masses, num_parts, size);
    //MPI_Barrier( MPI_COMM_WORLD);
    std::cout << "Init simulation ended." << std::endl;
    

    auto init_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff_1 = init_time - start_time;
    double seconds_1 = diff_1.count();
    std::cout << "initialization Time = " << seconds_1 << " seconds\n";

    */

   std::cout << rank << " -> displ_loc: " << displ_loc << std::endl;

    if(!rank) save(fsave, parts_pos, num_parts, size);
    //for nel tempo: non parallelizzare
    for (int step = 0; step < nsteps; ++step) {
        //MPI_Barrier( MPI_COMM_WORLD);
        
        simulate_one_step(parts_pos, parts_vel_acc_loc, masses, num_parts, num_loc, displ_loc, size, rank);
        
        MPI_Barrier( MPI_COMM_WORLD);
        MPI_Allgatherv( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL , &parts_pos[0] , &sizes[0] , &displs[0] , mpi_part_pos_type , MPI_COMM_WORLD);
        
        //MPI_Allgather( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL ,  &parts_pos[0] , num_loc, mpi_part_pos_type , MPI_COMM_WORLD);
        //MPI_Barrier( MPI_COMM_WORLD);
        
        // Save state if necessary
        if(rank==0)
        {
            //std::ofstream fsave(savename);
            if (fsave.good() && (step % savefreq) == 0) {
                //gather_for_save(parts, masses, num_parts, rank, size );
                //std::cout << "rank " << rank << " : saving at step " << step <<std::endl;
                save(fsave, parts_pos, num_parts, size);
            }
            
            if(step > 0){
                if (step%10 == 0){
                fflush(stdout);
                printf("[ %d% ]\r", (int)(step*100/nsteps));
                }
                
            }
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
