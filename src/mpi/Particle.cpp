#include "common.h"
#include "Particle.hpp"

#include <memory>
#include <stdexcept>
#include <cmath>


using namespace common_h;

    void particle_pos:: move(const double size, particle_vel_acc& parts_vel_acc ) {
        // Slightly simplified Velocity Verlet integration
        // Conserves energy better than explicit Euler method
        parts_vel_acc.vx += parts_vel_acc.ax * dt;
        parts_vel_acc.vy += parts_vel_acc.ay * dt;
        parts_vel_acc.vz += parts_vel_acc.az * dt;
        this->x += parts_vel_acc.vx * dt;
        this->y += parts_vel_acc.vy * dt;
        this->z += parts_vel_acc.vz * dt;

        // Bounce from walls
        while (this->x < 0 || this->x > size) {
            this->x = this->x < 0 ? -this->x : 2 * size - this->x;
            parts_vel_acc.vx = -parts_vel_acc.vx;
        }

        while (this->y < 0 || this->y > size) {
            this->y = this->y < 0 ? -this->y : 2 * size - this->y;
            parts_vel_acc.vy = -parts_vel_acc.vy;
        }

        while (this->z < 0 || this->z > size) {
            this->z = this->z < 0 ? -this->z : 2 * size - this->z;
            parts_vel_acc.vz = -parts_vel_acc.vz;
        }
    };


