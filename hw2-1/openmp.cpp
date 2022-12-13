#include "common.h"
#include <cmath>
#include <vector>
#include <random>

    // Apply the force from neighbor to particle
    void particle:: apply_force(particle& neighbor) {
        // Calculate Distance
        double dx = neighbor.x - this->x;
        double dy = neighbor.y - this->y;
        double r2 = dx * dx + dy * dy;

        // Check if the two particles should interact
        if (r2 > cutoff * cutoff)
            return;

        r2 = fmax(r2, min_r * min_r);
        double r = sqrt(r2);

        // Very simple short-range repulsive force
        double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
        this->ax += coef * dx;
        this->ay += coef * dy;
    }

    // Integrate the ODE

    void particle:: move(double size) {
        // Slightly simplified Velocity Verlet integration
        // Conserves energy better than explicit Euler method
        this->vx += this->ax * dt;
        this->vy += this->ay * dt;
        this->x += this->vx * dt;
        this->y += this->vy * dt;

        // Bounce from walls
        while (this->x < 0 || this->x > size) {
            this->x = this->x < 0 ? -this->x : 2 * size - this->x;
            this->vx = -this->vx;
        }

        while (this->y < 0 || this->y > size) {
            this->y = this->y < 0 ? -this->y : 2 * size - this->y;
            this->vy = -this->vy;
        }
    }



void simulate_one_step(std::vector<particle>& parts,int num_parts, double size) {
    // Compute Forces
    //int num_parts = parts.size();
	#pragma omp for schedule(dynamic)
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
        for (int j = 0; j < num_parts; ++j) {
            parts[i].apply_force(parts[j]);
        }
    }
	#pragma omp barrier

    // Move Particles
	#pragma omp for schedule(dynamic)
    for (int i = 0; i < num_parts; ++i) {
        parts[i].move(size);
    }
}
