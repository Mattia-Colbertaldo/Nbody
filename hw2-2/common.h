#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>
#include <mpi.h>

// Program Constants
#define nsteps   1000
#define savefreq 1
#define density  0.0005
#define cutoff   0.01
#define min_r    (cutoff / 100)
#define dt       0.0005
#define G        0.00000000000667

// Particle Data Structure: used in OPENMP
class particle_pos {
    public:

    particle_pos(){};
    particle_pos(const double x, const double y, const double z) : x(x), y(y), z(z){};
    double x;  // Position X
    double y;  // Position Y
    double z;
    
};


class particle_vel_acc {
    public:

    double vx; // Velocity X
    double vy; // Velocity Y
    double vz;
    double ax; // Acceleration X
    double ay; // Acceleration Y
    double az;
    particle_vel_acc(){}; 
    particle_vel_acc(const double vx, const double vy, const double vz) : vx(vx), vy(vy), ax(0.), ay(0.), az(0.){};
    void apply_force(particle_pos & me, particle_pos& neighbor, float mass);
    void move(particle_pos & me, double size);
};

#endif