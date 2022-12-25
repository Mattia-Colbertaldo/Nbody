#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>
#include <random>

// Program Constants
#define nsteps   1000
#define savefreq 1
#define density  0.0005
#define cutoff   0.01
#define min_r    (cutoff / 100)
#define dt       0.0005
#define G        0.00000000000667

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
        
        void apply_force(particle& neighbor);
        void move(double size);
        

};

#endif