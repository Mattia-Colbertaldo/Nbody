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





/******         SAVING       *******/

void save(std::ofstream& fsave, std::vector<particle_pos>& parts, int num_parts, double size) {

    if (first) {
        fsave << num_parts << " " << size << "\n";
        first = false;
    }
    

    for (int i = 0; i < num_parts; ++i) {
        fsave << parts[i].x << " " << parts[i].y << " " << parts[i].z <<"\n";
    }

    fsave << std::endl;
}





/******         PARTICLE INITIALIZATION       *******/


void init_particles(std::vector<particle_vel_acc>& parts_vel_acc_loc,
                    std::vector<float>& masses,
                    int num_parts, double size,int part_seed, int rank,
                    std::vector<particle_pos>&parts_pos, int num_loc, 
                     std::vector<int>& displs, std::vector<int>& sizes) {

    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());
    
    int sx = (int)ceil(sqrt((double)num_parts));
    int sy = (num_parts + sx - 1) / sx;
    int sz = sy;


    ////// MPI STRUCT /////////////

    const int nitems=6;
    int          blocklengths[6] = {1,1,1,1,1,1};
    MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    //MPI_Datatype mpi_part_vel_acc_type;
    MPI_Aint     offsets[6];

    offsets[0] = offsetof(particle_vel_acc, vx);
    offsets[1] = offsetof(particle_vel_acc, vy);
    offsets[2] = offsetof(particle_vel_acc, vz);
    offsets[3] = offsetof(particle_vel_acc, ax);
    offsets[4] = offsetof(particle_vel_acc, ay);
    offsets[5] = offsetof(particle_vel_acc, az);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_part_vel_acc_type);
    MPI_Type_commit(&mpi_part_vel_acc_type);
    
    ///////////////////////////////

    ////// MPI STRUCT /////////////
    
    const int nitems1=3;
    int          blocklengths1[3] = {1,1,1};
    MPI_Datatype types1[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    //MPI_Datatype mpi_part_pos_type;
    MPI_Aint     offsets1[3];

    offsets1[0] = offsetof(particle_pos, x);
    offsets1[1] = offsetof(particle_pos, y);
    offsets1[2] = offsetof(particle_pos, z);

    MPI_Type_create_struct(nitems1, blocklengths1, offsets1, types1, &mpi_part_pos_type);
    MPI_Type_commit(&mpi_part_pos_type);
    
    ///////////////////////////////


    std::vector<int> shuffle(num_parts);
    for (int i = 0; i < shuffle.size(); ++i) {
        shuffle[i] = i;
    }
    
    std::vector<particle_vel_acc> parts_vel_acc_temp(num_parts);

    if(rank==0){


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
            parts_pos[i].z = parts_pos[i].x * parts_pos[i].y;
            
            // Assign random velocities within a bound
            std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
            parts_vel_acc_temp[i].vx = rand_real(gen);
            parts_vel_acc_temp[i].vy = rand_real(gen);
            parts_vel_acc_temp[i].vz = rand_real(gen);

            std::uniform_real_distribution<float> rand_mass(0.001, 0.1);
            float m = rand_mass(gen);
            masses[i]=m;

        }

        //Saving my data (i'm rank 0)
        for(int i=0; i<num_loc; i++){
            parts_vel_acc_loc[i].vx = parts_vel_acc_temp[i].vx;
            parts_vel_acc_loc[i].vy = parts_vel_acc_temp[i].vy;
            parts_vel_acc_loc[i].vz = parts_vel_acc_temp[i].vz;
        }
         
    }
    MPI_Scatterv( &parts_vel_acc_temp[0] , &sizes[0] , &displs[0], mpi_part_vel_acc_type ,
                  &parts_vel_acc_loc[0] , sizes[rank] , mpi_part_vel_acc_type , 0, MPI_COMM_WORLD);
    MPI_Bcast( &masses[0] , num_parts , MPI_FLOAT , 0 , MPI_COMM_WORLD);
    MPI_Bcast(&parts_pos[0] , num_parts, mpi_part_pos_type , 0 , MPI_COMM_WORLD);
}







/******         ONE STEP SIMULATION       *******/


void simulate_one_step(std::vector<particle_pos>& parts_pos, std::vector<particle_vel_acc>& parts_vel_acc_loc,
                       std::vector<float>& masses, int num_parts, int num_loc, int displ_loc,  double size, int rank){

    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `mpi_rank`
    // in this way we split the load in the most equilibrate way
    // Ogni processore aggiorna le particelle nel range [mpi_rank*N, (mpi_rank+1)*N).
    // Notate che per utilizzare apply_force e move vi servono posizione, velocità e massa
    // delle particelle in [mpi_rank*N, (mpi_rank+1)*N) e solo posizione e massa delle particelle in [0, N)

    for (int i = 0; i < num_loc; ++i) {
        parts_vel_acc_loc[i].ax = parts_vel_acc_loc[i].ay = 0;
        for (int j = 0; j < num_parts; ++j) {
            //OK;
            parts_vel_acc_loc[i].apply_force(parts_pos[i+displ_loc], parts_pos[j], masses[j]);
            //OK;
        }
    }

    // Move Particles
	
    for (int i = 0; i < num_loc; ++i) {
        //OK;
        parts_vel_acc_loc[i].move(parts_pos[i+displ_loc] , size);
        //OK;
        
    }
     
    
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
            //std::cout << "Size of process " << i << " = " << sizes[i] << std::endl;
            displs[i+1]= displs[i]+sizes[i];
            //std::cout << "Displ of process " << i << " = " << displs[i] << std::endl;
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
    num_loc=sizes[rank];
    displ_loc=displs[rank];
    std::vector<particle_vel_acc> parts_vel_acc_loc(num_loc);
    std::vector<particle_pos> parts_pos(num_parts);
    
    std::vector<float> masses(num_parts);
    
    if (!rank) std::cout << "Trying to init particles..." << std::endl;
    init_particles(parts_vel_acc_loc, masses, num_parts, size, part_seed, rank, parts_pos, num_loc, displs, sizes); 
    
    
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

   if(!rank){
        auto init_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff_1 = init_time - start_time;
        double seconds_1 = diff_1.count();
        std::cout << "initialization Time = " << seconds_1 << " seconds\n";
    }


    if(!rank) save(fsave, parts_pos, num_parts, size);
    //for nel tempo: non parallelizzare
    for (int step = 0; step < nsteps; ++step) {
        
        simulate_one_step(parts_pos, parts_vel_acc_loc, masses, num_parts, num_loc, displ_loc, size, rank);
        
        MPI_Barrier( MPI_COMM_WORLD);

        // Allgather delle posizioni, in questo modo aggiorno la posizione di tutte le particelle per tutti i processori.
        // Non serve comunicare velocità e accelerazione visto che sono necessarie solo localmente.
        MPI_Allgatherv( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL , &parts_pos[0] , &sizes[0] , &displs[0] , mpi_part_pos_type , MPI_COMM_WORLD);
        
        // Save state if necessary
        if(rank==0)
        {
            if (fsave.good() && (step % savefreq) == 0) {
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
