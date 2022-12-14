#include "common.h"
#include <cmath>
#include <mpi.h>
#include <iostream>

int mpi_rank;


#define OK std::cout << "At mpi:" << __LINE__ << " from process " << mpi_rank << std::endl


// Apply the force from neighbor to particle
void particle_vel_acc :: apply_force(const particle_pos& me, const particle_pos& neighbor, float mass_neigh) {
    // Calculate Distance
    double dx = neighbor.x - me.x;
    double dy = neighbor.y - me.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;
    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = mass_neigh* (1 - cutoff / r) / r2;
    this->ax += coef * dx;
    this->ay += coef * dy;
}

// Integrate the ODE
void particle_vel_acc:: move(particle_pos& pos ,double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    this->vx += this->ax * dt;
    this->vy += this->ay * dt;
    pos.x += this->vx * dt;
    pos.y += this->vy * dt;
    
    // Bounce from walls
    while (pos.x < 0 || pos.x > size) {
        pos.x = pos.x < 0 ? -pos.x : 2 * size - pos.x;
        this->vx = -this->vx;
    }
    
    while (pos.y < 0 || pos.y > size) {
        pos.y = pos.y < 0 ? -pos.y : 2 * size - pos.y;
        this->vy = -this->vy;
    }
}


void init_simulation(std::vector<particle_pos>& parts,std::vector<float>& masses,int num_parts, double size) {

	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here

    // broadcast the global size to prepare the local buffers
    
    

    // since each chunk has a different size we have to prepare buffers with
    // sizes and displacements of the chunk we have to send

    

}

void simulate_one_step(std::vector<particle_pos>& parts_pos, std::vector<particle_vel_acc>& parts_vel_acc_loc, std::vector<float>& masses, int num_parts, int num_loc, int displ_loc,  double size, int rank){
    mpi_rank = rank;

    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `mpi_rank`
    // in this way we split the load in the most equilibrate way

   //Ogni processore aggiorna le particelle nel range [mpi_rank*N, (mpi_rank+1)*N). Notate che per utilizzare apply_force e move vi servono posizione, velocità e massa delle particelle in [mpi_rank*N, (mpi_rank+1)*N) e solo posizione e massa delle particelle in [0, N)
    for (int i = 0; i < num_loc; ++i) {
        parts_vel_acc_loc[i].ax = parts_vel_acc_loc[i].ay = 0;
        for (int j = 0; j < num_parts; ++j) {
            //OK;
            parts_vel_acc_loc[i].apply_force(parts_pos[i+displ_loc], parts_pos[j], masses[j]);
            //OK;
        }
    }

    //MPI_Barrier( MPI_COMM_WORLD );

    
    


    // Move Particles
	
    for (int i = 0; i < num_loc; ++i) {
        //OK;
        parts_vel_acc_loc[i].move(parts_pos[i+displ_loc] , size);
        //OK;
        
    }
    // Allgather delle posizioni, in questo modo aggiorno la posizione di tutte le particelle per tutti i processori. Non serve comunicare velocità e accelerazione visto che sono necessarie solo localmente. 
    
}


/*
void gather_for_save(particle_t* parts, ma int num_parts, double size, int mpi_rank, int num_procs) {
    // Write this function such that at the end of it, the master (mpi_rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.
}
*/
