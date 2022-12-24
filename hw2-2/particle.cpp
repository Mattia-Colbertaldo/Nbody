#include "common.h"
#include <cmath>
#include <mpi.h>
#include <iostream>
int mpi_rank;


#define OK std::cout << "At mpi:" << __LINE__ << " from process " << mpi_rank << std::endl


// Apply the force from neighbor to particle
void particle_vel_acc :: apply_force(particle_pos& me, particle_pos& neighbor, float mass_neigh) {
    // Calculate Distance
    double dx = neighbor.x - me.x;
    double dy = neighbor.y - me.y;
    double dz = neighbor.z - me.z;
    double r2 = dx * dx + dy * dy + dz * dz;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;
    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = mass_neigh* (1 - cutoff / r) / r2;
    this->ax += coef * dx;
    this->ay += coef * dy;
    this->az += coef * dz; 
}

// Integrate the ODE
void particle_vel_acc:: move(particle_pos& pos ,double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    this->vx += this->ax * dt;
    this->vy += this->ay * dt;
    this->vz += this->az * dt;
    pos.x += this->vx * dt;
    pos.y += this->vy * dt;
    pos.z += this->vz * dt;
    
    // Bounce from walls
    while (pos.x < 0 || pos.x > size) {
        pos.x = pos.x < 0 ? -pos.x : 2 * size - pos.x;
        this->vx = -this->vx;
    }
    
    while (pos.y < 0 || pos.y > size) {
        pos.y = pos.y < 0 ? -pos.y : 2 * size - pos.y;
        this->vy = -this->vy;
    }

    while (pos.z < 0 || pos.z > size) {
        pos.z = pos.z < 0 ? -pos.z : 2 * size - pos.z;
        this->vz = -this->vz;
    }
}