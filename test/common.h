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
typedef struct particle_pos {
    double x;  // Position X
    double y;  // Position Y
} particle_pos;

typedef struct particle_vel_acc {
    double vx; // Velocity X
    double vy; // Velocity Y
    double ax; // Acceleration X
    double ay; // Acceleration Y
} particle_vel_acc;

typedef struct particle_pos_vel {
    double x;  // Position X
    double y;  // Position Y
    double vx; // Velocity X
    double vy; // Velocity Y
} particle_pos_vel;




//for MPI
void init_simulation(std::vector<particle_pos>& parts,std::vector<float>& masses,int num_parts, double size);
void simulate_one_step(std::vector<particle_pos_vel>& parts_pos_vel_loc, std::vector<particle_pos>& parts_pos, std::vector<particle_vel_acc>& parts_vel_acc_loc, std::vector<float>& masses, int num_parts, double size, std::vector<int> &sizes,std::vector<int> & displs);
/*void gather_for_save(std::vector<particle_mpi> parts, std::vector<float>& masses, int num_parts, int rank, double size);*/

#endif