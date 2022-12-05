#include "common.h"
#include <cmath>
#include <mpi.h>

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
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
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}


void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here

    // broadcast the global size to prepare the local buffers
    
    MPI_Bcast(&num_parts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // the local size is `n / size` plus 1 if the reminder `n % size` is greater than `rank`
    // in this way we split the load in the most equilibrate way
    const auto local_size = n / size + (n % size > rank);
    std::vector<double> local_x(local_size), local_y(local_size), local_vx(local_size), local_vy(local_size), local_ax(local_size), local_ay(local_size);

    // since each chunk has a different size we have to prepare buffers with
    // sizes and displacements of the chunk we have to send
    std::vector<int> sizes(size), displs(size + 1);
    for (int i = 0; i < size; ++i) {
        sizes[i] = n / size + (n % size > i);
        displs[i + 1] = displs[i] + sizes[i];
    }

    // scatter the chunks
    
    MPI_Bcast( parts , num_parts , MPI_DOUBLE , 0 , MPI_COMMON_WORLD); //FLAG
    
    MPI_Scatterv(parts->vx, sizes, displs, MPI_DOUBLE,
                local_vx, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(parts->vy, sizes, displs, MPI_DOUBLE,
                local_vy, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(parts->ax, sizes, displs, MPI_DOUBLE,
                local_ax, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(parts->ay, sizes, displs, MPI_DOUBLE,
                local_ay, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces
    int size_t;
    MPI_Comm_size(MPI_COMM_WORLD, &size_t);
	for(int t=0; t<size; i++){
        for (int i = displs[t]; i < sizes[t]; ++i) {
            parts[i].ax = parts[i].ay = 0;
            for (int j = 0; j < num_parts; ++j) {
                apply_force(parts[i], parts[j]);
            }
        }
    }
    
	

    // Move Particles
	
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }

    MPI_Allgather( parts , size , MPI_ , parts ,size , MPI_ , MPI_COMM_WORLD); //FLAG
    
}
