#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>
#include <mpi.h>

// Program Constants
#define nsteps   1000
#define savefreq 10
#define density  0.0005
#define cutoff   0.01
#define min_r    (cutoff / 100)
#define dt       0.0005
#define G        0.00000000000667

// Particle Data Structure: used in OPENMP
class particle_pos {
    public:
    particle_pos(){};
    particle_pos(const double x, const double y){
    this->x = x;
    this->y = y;
    };
        
    double x;  // Position X
    double y;  // Position Y
    
};


class particle_vel_acc {
    public:
    
    double vx; // Velocity X
    double vy; // Velocity Y
    double ax; // Acceleration X
    double ay; // Acceleration Y
    particle_vel_acc(){};
            
    particle_vel_acc(const double vx, const double vy){
    this->vx = vx;
    this->vy = vy;
    this->ax = 0.;
    this->ay = 0.;
    };
    void apply_force(particle_pos & me, particle_pos& neighbor, float mass);
    void move(particle_pos & me, double size);
};



//for MPI
void init_simulation(std::vector<particle_pos>& parts,std::vector<float>& masses,int num_parts, double size);
void simulate_one_step(std::vector<particle_pos>& parts_pos, std::vector<particle_vel_acc>& parts_vel_acc_loc, std::vector<float>& masses, int num_parts, int num_loc, int displ_loc, double size, int rank);

/*void gather_for_save(std::vector<particle_mpi> parts, std::vector<float>& masses, int num_parts, int rank, double size);*/

#endif