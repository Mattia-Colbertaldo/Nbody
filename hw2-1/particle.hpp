#ifndef HH_PARTICLE_HH
#define HH_PARTICLE_HH
#include <memory>
#include <vector>
#include <random>
#include "Force.hpp"


// Particle Data Structure: used in OPENMP

class particle {
    public:
        
        double x;  // Position X
        double y;  // Position Y
        double z;  // Position Z
        double vx; // Velocity X
        double vy; // Velocity Y
        double vz; // Velocity Z
        double ax; // Acceleration X
        double ay; // Acceleration Y
        double az; // Acceleration Z
        float mass;
        particle(){};
        
        particle(const double x, const double y, const double z,
                const double vx, const double vy, const double vz,
                const float m) : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz),
                ax(0.), ay(0.), az(0.), mass(m){};
        
        void apply_force(const particle& neighbor, const std::shared_ptr<AbstractForce> force);
        void move(double size);
        

};



#endif