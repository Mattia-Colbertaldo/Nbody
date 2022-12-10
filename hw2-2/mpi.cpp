#include "common.h"
#include <cmath>
#include <mpi.h>
#include <iostream>

int mpi_rank;

#define OK std::cout << "At mpi:" << __LINE__ << " from process " << mpi_rank << std::endl


// Apply the force from neighbor to particle
void apply_force(particle_vel_acc& particle_vel_acc_loc,  particle_pos& particle,  particle_pos& neighbor, float mass_neigh) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = mass_neigh* (1 - cutoff / r) / r2;
    particle_vel_acc_loc.ax += coef * dx;
    particle_vel_acc_loc.ay += coef * dy;
}

// Integrate the ODE
void move(particle_vel_acc& particle_vel_acc_loc, particle_pos& pos ,double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    particle_vel_acc_loc.vx += particle_vel_acc_loc.ax * dt;
    particle_vel_acc_loc.vy += particle_vel_acc_loc.ay * dt;
    pos.x += particle_vel_acc_loc.vx * dt;
    pos.y += particle_vel_acc_loc.vy * dt;

    // Bounce from walls
    while (pos.x < 0 || pos.x > size) {
        pos.x = pos.x < 0 ? -pos.x : 2 * size - pos.x;
        particle_vel_acc_loc.vx = -particle_vel_acc_loc.vx;
    }

    while (pos.y < 0 || pos.y > size) {
        pos.y = pos.y < 0 ? -pos.y : 2 * size - pos.y;
        particle_vel_acc_loc.vy = -particle_vel_acc_loc.vy;
    }
}


void init_simulation(std::vector<particle_pos>& parts,std::vector<float>& masses,int num_parts, double size) {
    //int num_parts = parts.size();

	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here

    // broadcast the global size to prepare the local buffers
    
    

    // since each chunk has a different size we have to prepare buffers with
    // sizes and displacements of the chunk we have to send

    

}

void simulate_one_step( std::vector<particle_pos>& parts_pos, std::vector<particle_vel_acc>& parts_vel_acc_loc, std::vector<float>& masses, int num_parts, int num_loc, double size) {
   
    // Compute Forces
    //int num_parts = parts.size();
    int mpi_size;
    MPI_Comm_size( MPI_COMM_WORLD , &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    //std::cout << "I'm process " << mpi_rank << std::endl;
   
    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `mpi_rank`
    // in this way we split the load in the most equilibrate way


    /*
    for(int t=0; t<num_parts; t++){
        std::cout << t << "-" << parts[t].x << " " << parts[t].y << " ";
    }
    std::cout << std::endl;
    */
    /*
    for(int t=0; t<sizes[mpi_rank]/2; t++){
        std::cout << t << "-" << loc_parts[t].x << " " << loc_parts[t].y << " ";
    }
    std::cout << "From process " << mpi_rank << std::endl;
    */

   //Ogni processore aggiorna le particelle nel range [mpi_rank*N, (mpi_rank+1)*N). Notate che per utilizzare apply_force e move vi servono posizione, velocità e massa delle particelle in [mpi_rank*N, (mpi_rank+1)*N) e solo posizione e massa delle particelle in [0, N)
    for (int i = 0; i < num_loc; ++i) {

        parts_vel_acc_loc[i].ax = parts_vel_acc_loc[i].ay = 0;
        
        for (int j = 0; j < num_parts; ++j) {
            apply_force(parts_vel_acc_loc[i], parts_pos[i], parts_pos[j], masses[j]);
        }
        
    }



    // Move Particles
	
    for (int i = 0; i < num_loc; ++i) {
        move(parts_vel_acc_loc[i],parts_pos[i+num_loc*mpi_rank] , size);
    }
    // Allgather delle posizioni, in questo modo aggiorno la posizione di tutte le particelle per tutti i processori. Non serve comunicare velocità e accelerazione visto che sono necessarie solo localmente. 
    
    MPI_Allgather( MPI_IN_PLACE , 0 , MPI_DATATYPE_NULL ,  &parts_pos[0] , num_loc*2 , MPI_DOUBLE , MPI_COMM_WORLD);

}


/*
void gather_for_save(particle_t* parts, ma int num_parts, double size, int mpi_rank, int num_procs) {
    // Write this function such that at the end of it, the master (mpi_rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.
}
*/