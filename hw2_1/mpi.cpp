#include "common.h"
#include <cmath>
#include <mpi.h>

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor, float mass_neigh) {
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
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p,  std::vector<particle_mpi>& loc_parts, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    loc_parts.x += p.vx * dt;
    loc_parts.y += p.vy * dt;

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


void init_simulation(std::vector<particle_mpi>& parts,int num_parts, double size) {
    //int num_parts = parts.size();

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

    // todo: optimize this by broadcasting directly the particle
    std::vector<double> x(num_parts), y(num_parts), vx(num_parts), vy(num_parts), ax(num_parts), ay(num_parts);
    for(int i=0; i<num_parts; i++){
        x.emplace_back(parts[i].x);
        y.emplace_back(parts[i].y);
        vx.emplace_back(parts[i].vx);
        vy.emplace_back(parts[i].vy);
        ax.emplace_back(parts[i].ax);
        ay.emplace_back(parts[i].ay);
    }

    // scatter the chunks
    
    auto partsdata = parts.data()

    std::vector<particle_mpi> loc_parts(size[rank]);

    MPI_Bcast(masses , size[rank] , MPI_FLOAT , 0 , MPI_COMM_WORLD); //FLAG BCAST MASSES
    
   
    MPI_Scatterv(parts, size[rank]*2 , displs, MPI_DOUBLE,
                loc_parts, size[rank]*2 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

}

void simulate_one_step( std::vector<particle_mpi>& loc_parts,int num_parts, double size) {
    // Compute Forces
    //int num_parts = parts.size();
    int size_t;
    MPI_Comm_size(MPI_COMM_WORLD, &size_t);
	for(int t=0; t<size; i++){
        for (int i = displs[t]; i < sizes[t]; ++i) {
            std::vector<particle_t>particle_t part_acc;
            part_acc[i].x=loc_parts[i].x;
            part_acc[i].y=loc_parts[i].y; 
            part_acc[i].ax = part_acc[i].ay = 0;
            for (int j = 0; j < num_parts; ++j) {
                apply_force(part_acc[i], part_acc[j], masses[j]);
            }
        }
    }
    
	

    // Move Particles
	
    for (int i = 0; i < num_parts; ++i) {
        move(part_acc[i], loc_parts[i] , size);
    }

   
    MPI_Allgather( loc_parts , size*2 , MPI_DOUBLE , parts , size*2 , MPI_DOUBLE , MPI_COMMM_WORLD); //FLAG
    
}
