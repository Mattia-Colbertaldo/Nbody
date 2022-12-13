#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>

// Program Constants
#define nsteps   1000
#define savefreq 10
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
        double vx; // Velocity X
        double vy; // Velocity Y
        double ax; // Acceleration X
        double ay; // Acceleration Y
        float mass;
        void apply_force(particle& neighbor);
        void move(double size);

};


// Simulation routine
void simulate_one_step(std::vector<particle>& parts, int num_parts, double size);


#endif