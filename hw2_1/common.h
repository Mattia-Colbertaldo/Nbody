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
typedef struct particle_t {
    double x;  // Position X
    double y;  // Position Y
    double vx; // Velocity X
    double vy; // Velocity Y
    double ax; // Acceleration X
    double ay; // Acceleration Y
} particle_t;

// Particle Data Structure: used in MPI
typedef struct particle_mpi {
    
    double x;    // Position X
    double y;    // Position Y
} particle_mpi;

//for MPI
void init_simulation(std::vector<particle_mpi>& parts, std::vector<float>& masses, int num_parts, double size);
void simulate_one_step(std::vector<particle_mpi>& parts, int num_parts, double size);

// Simulation routine
void init_simulation(std::vector<particle_t>& parts, std::vector<float>& masses, int num_parts,double size);
void simulate_one_step(std::vector<particle_t>& parts, int num_parts, double size);


#endif