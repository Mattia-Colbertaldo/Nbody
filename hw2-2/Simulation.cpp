#include "Particle.hpp"
#include "Simulation.hpp"
#include <memory>
#include <stdexcept>
#include <cmath>
#include <random>
#include <iostream>


/******         PARTICLE INITIALIZATION       *******/


void Simulation :: init_particles(
                    const int num_parts, const double size,const int part_seed,
                    const int num_loc, 
                    const std::vector<int>& displs, const std::vector<int>& sizes) {

    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());
    std::mt19937 gen2(part_seed ? part_seed : rd());
    
    int sx = (int)ceil(std::sqrt((double)num_parts));
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

    offsets1[0] = offsetof(this->part_pos, x);
    offsets1[1] = offsetof(this->part_pos, y);
    offsets1[2] = offsetof(this->part_pos, z);

    MPI_Type_create_struct(nitems1, blocklengths1, offsets1, types1, &mpi_part_pos_type);
    MPI_Type_commit(&mpi_part_pos_type);
    
    ///////////////////////////////


    std::vector<int> shuffle(num_parts);
    for (int i = 0; i < shuffle.size(); ++i) {
        shuffle[i] = i;
    }
    
    std::vector<part_vel_acc> parts_vel_acc_temp(num_parts);

    if(rank==0){


        // initialize local vector of positions and velocities (this->part_pos_vel_loc)
        // also fill the local part of vector of positions (this->part_pos) and then allgather it
        // --> result : all have local values of positions and velocities in this->part_pos_vel_loc and the positions of ALL particles in parts.this->part_pos
        for (int i = 0; i < num_parts; ++i) {
            // Make sure particles are not spatially sorted
            std::uniform_int_distribution<int> rand_int(0, num_parts - i - 1);
            int j = rand_int(gen);
            int k = shuffle[j];
            shuffle[j] = shuffle[num_parts - i - 1];
            // Distribute particles evenly to ensure proper spacing
            
            this->part_pos[i].x = size * (1. + (k % sx)) / (1 + sx);
            this->part_pos[i].y = size * (1. + (k / sx)) / (1 + sy);
            this->part_pos[i].z = this->part_pos[i].x * this->part_pos[i].y;
            
            // Assign random velocities within a bound
            std::uniform_real_distribution<double> rand_real(-1.0, 1.0);
            parts_vel_acc_temp[i].vx = rand_real(gen);
            parts_vel_acc_temp[i].vy = rand_real(gen);
            parts_vel_acc_temp[i].vz = rand_real(gen);

            std::uniform_real_distribution<double> rand_mass(0.001, 0.1);
            double m = rand_mass(gen);
            this->masses[i]=m;

            std::uniform_real_distribution<double> rand_charge(-1.0, 1.0);
            double charge=rand_charge(gen2) * 1e-19;
            this->charges[i]=charge;

        }

        //Saving my data (i'm rank 0)
        for(int i=0; i<num_loc; i++){
            parts_vel_acc_loc[i].vx = parts_vel_acc_temp[i].vx;
            parts_vel_acc_loc[i].vy = parts_vel_acc_temp[i].vy;
            parts_vel_acc_loc[i].vz = parts_vel_acc_temp[i].vz;
        }
         
    }
    MPI_Scatterv( &parts_vel_acc_temp[0] , &sizes[0] , &displs[0], mpi_part_vel_acc_type ,
                  &(this->parts_vel_acc_loc[0]) , sizes[rank] , mpi_part_vel_acc_type , 0, MPI_COMM_WORLD);
    MPI_Bcast( &this->masses[0] , num_parts , MPI_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &this->charges[0] , num_parts , MPI_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Bcast(&this->part_pos[0] , num_parts, mpi_part_pos_type , 0 , MPI_COMM_WORLD);
};


    



/******         ONE STEP SIMULATION       *******/


void Simulation :: simulate_one_step(int num_parts, int num_loc, int displ_loc,  double size, int rank, const std::shared_ptr<AbstractForce>& force ){

    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `mpi_rank`
    // in this way we split the load in the most equilibrate way
    // Ogni processore aggiorna le particelle nel range [mpi_rank*N, (mpi_rank+1)*N).
    // Notate che per utilizzare apply_force e move vi servono posizione, velocità e massa
    // delle particelle in [mpi_rank*N, (mpi_rank+1)*N) e solo posizione e massa delle particelle in [0, N)

    for (int i = 0; i < num_loc; ++i) {
        this->parts_vel_acc_loc[i].ax = this->parts_vel_acc_loc[i].ay = this->parts_vel_acc_loc[i].az= 0.;
        for (int j = 0; j < num_parts; ++j) {
            if(i+displ_loc != j) force.force_application(this->parts_pos, this->parts_vel_acc_loc, this->masses[j], this->charges[i], this->charges[j], i, j);
        }
    }

    // Move Particles
    for (int i = 0; i < num_loc; ++i) {
        this->parts_vel_acc_loc[i].move(this->parts_pos[i+displ_loc] , size);
        
    }
     
};