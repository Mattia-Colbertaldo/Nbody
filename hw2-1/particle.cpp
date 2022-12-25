
#include "common.h"

#include <memory>
#include <cmath>

    // Apply the force from neighbor to particle
    void particle:: apply_force(const particle& neighbor, const std::string forcename) {
        const std::unordered_map<std::string, std::shared_ptr<AbstractForce>> fmap =
        {
            {"repulsive", std::make_shared<RepulsiveForce>() },
            {"gravitational", std::make_shared<GravitationalForce>() },
            {"coulomb", std::make_shared<CoulombForce>() },
        };
        
        std::shared_ptr<AbstractForce> force; 
        try {
            force = fmap.at(forcename);      // vector::at throws an out-of-range
        }
        catch (const std::out_of_range& oor) {
            std::cerr << "Default Force chosen "<< '\n';
            force =  fmap.at("repulsive");
        }
        force->force_application (*this, neighbor);
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



