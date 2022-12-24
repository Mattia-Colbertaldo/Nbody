#include "common.h"
#include <cmath>
#include <iostream>

        
        void particle:: apply_force(particle& neighbor) {
        // Calculate Distance
        double dx = neighbor.x - this->x;
        double dy = neighbor.y - this->y;
        double dz = neighbor.z - this->z;
        double r2 = dx * dx + dy * dy + dz* dz;

        // Check if the two particles should interact
        if (r2 > cutoff * cutoff)
            return;

        r2 = fmax(r2, min_r * min_r);
        double r = sqrt(r2);

        // Very simple short-range repulsive force
        double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
        this->ax += coef * dx;
        this->ay += coef * dy;
        this->az += coef * dz;
    }

    // Integrate the ODE

    void particle:: move(double size) {
        // Slightly simplified Velocity Verlet integration
        // Conserves energy better than explicit Euler method
        this->vx += this->ax * dt;
        this->vy += this->ay * dt;
        this->vz += this->az * dt;
        this->x += this->vx * dt;
        this->y += this->vy * dt;
        this->z += this->vz * dt;

        // Bounce from walls
        while (this->x < 0 || this->x > size) {
            this->x = this->x < 0 ? -this->x : 2 * size - this->x;
            this->vx = -this->vx;
        }

        while (this->y < 0 || this->y > size) {
            this->y = this->y < 0 ? -this->y : 2 * size - this->y;
            this->vy = -this->vy;
        }

        while (this->z < 0 || this->z > size) {
            this->z = this->z < 0 ? -this->z : 2 * size - this->z;
            this->vz = -this->vz;
        }
    }


        // Apply the force from neighbor to particle
    

    // Integrate the ODE




void simulate_one_step(std::vector<particle>& parts, int num_parts, double size) {
    //int num_parts = parts.size();
    // Compute Forces
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = parts[i].az = 0.;
        for (int j = 0; j < num_parts; ++j) {
            parts[i].apply_force(parts[j]);
        }
    }


    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        parts[i].move(size);
    }
}


