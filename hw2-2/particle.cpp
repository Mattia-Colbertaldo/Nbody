#include "common.h"
#include "Particle.hpp"

#include <memory>
#include <stdexcept>
#include <cmath>


    // Integrate the ODE

    void Particle:: move(const double size) {
        // Slightly simplified Velocity Verlet integration
        // Conserves energy better than explicit Euler method
        this->part_vel_acc.vx += this->part_vel_acc.ax * dt;
        this->part_vel_acc.vy += this->part_vel_acc.ay * dt;
        this->part_vel_acc.vz += this->part_vel_acc.az * dt;
        this->part_pos.x += this->part_vel_acc.vx * dt;
        this->part_pos.y += this->part_vel_acc.vy * dt;
        this->part_pos.z += this->part_vel_acc.vz * dt;

        // Bounce from walls
        while (this->part_pos.x < 0 || this->part_pos.x > size) {
            this->part_pos.x = this->part_pos.x < 0 ? -this->part_pos.x : 2 * size - this->part_pos.x;
            this->part_vel_acc.vx = -this->part_vel_acc.vx;
        }

        while (this->part_pos.y < 0 || this->part_pos.y > size) {
            this->part_pos.y = this->part_pos.y < 0 ? -this->part_pos.y : 2 * size - this->part_pos.y;
            this->part_vel_acc.vy = -this->part_vel_acc.vy;
        }

        while (this->part_pos.z < 0 || this->part_pos.z > size) {
            this->part_pos.z = this->part_pos.z < 0 ? -this->part_pos.z : 2 * size - this->part_pos.z;
            this->part_vel_acc.vz = -this->part_vel_acc.vz;
        }
    };


