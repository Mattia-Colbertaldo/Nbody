#ifndef HH__SIMULATION__HH
#define HH__SIMULATION__HH

#include "PhysicalForce.hpp"
#include "Particle.hpp"
#include <vector>
#include <memory>
#include <mpi.h>

struct Simulation {
        public:

        Simulation(const int num_parts, const int num_loc, const int collision): collision(collision){
                this->parts_pos.resize(num_parts);
                this->parts_vel_acc_loc.resize(num_loc);
                this->masses.resize(num_parts);
                this->charges.resize(num_parts);
        };
        
        void simulate_one_step(int num_parts, int num_loc, int displ_loc, double size, int rank, const std::unique_ptr<AbstractForce>& force );

        MPI_Datatype init_particles(const int num_parts, const double size,const int part_seed,
                    const int num_loc, 
                    const std::vector<int>& displs, const std::vector<int>& sizes, const int rank);
        
        // private:
        std::vector<particle_pos> parts_pos;
        std::vector<particle_vel_acc> parts_vel_acc_loc;
        std::vector<double> masses;
        std::vector<double> charges;
        int collision;


};

#endif