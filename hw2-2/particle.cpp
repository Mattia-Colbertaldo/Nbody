#include "common.h"
#include <cmath>
#include <mpi.h>
#include <iostream>
int mpi_rank;


#define OK std::cout << "At mpi:" << __LINE__ << " from process " << mpi_rank << std::endl



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