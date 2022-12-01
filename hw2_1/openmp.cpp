#include "common.h"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

typedef struct force {
    double x;  // Force on X
    double y;  // Force on Y
} force;

std::vector<force> forces;
std::vector<std::vector<force>> loc_forces;

// Put any static global variables here that you will use throughout the simulation.

force apply_force(particle_t particle, particle_t neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff){
		force zero;
		zero.x = 0;
		zero.y = 0;
		return zero;
	}

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;

	force force_qk;
	force_qk.x = G * mass * mass / r2 * dx;
	force_qk.y = G * mass * mass / r2 * dy;
	return force_qk;

}

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here

	force zero;
	zero.x = 0;
	zero.y = 0;
	
	//forces[num_parts];
	//Vector of force.x and force.y for each particle
	forces = std::vector<force>(num_parts, zero);


	//loc_forces[omp_get_num_threads()][num_parts];
	loc_forces = std::vector<std::vector<force>>(omp_get_num_threads(), forces);

}

void simulate_one_step(particle_t* parts, int num_parts, double size) {

	//std::cout<<"Simulating one step"<<std::endl;

    // Write this function
	// Compute Forces
    #pragma omp parallel for schedule(dynamic) shared(num_parts, loc_forces)
	for (int q = 0; q<num_parts; ++q)
	{		
		#pragma omp parallel for shared(loc_forces)
		for(int k=q+1; k<num_parts; k++)
		{
			force force_qk;
			force_qk.x = 0;
			force_qk.y = 0;

			force_qk = apply_force(parts[q], parts[k]);

			int rank = omp_get_thread_num();
			
			loc_forces [ rank ] [ q ].x += force_qk.x ;
			loc_forces [ rank ] [ q ].y += force_qk.y ;

			loc_forces [ rank ] [ k ].x -= force_qk.x ;
			loc_forces [ rank ] [ k ].y -= force_qk.y ;
		}
	}
	#pragma omp parallel for schedule(dynamic)
	for (int q = 0; q<num_parts; ++q)
	{
		forces[q].x = 0;
		forces[q].y = 0;
		for(int t=0; t<omp_get_num_threads(); t++)
		{
			forces[q].x += loc_forces[t][q].x;
			forces[q].y += loc_forces[t][q].y;
		}
	}


}
