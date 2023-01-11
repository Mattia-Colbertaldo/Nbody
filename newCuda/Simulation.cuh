#ifndef HH__SIMULATION__HH
#define HH__SIMULATION__HH

#include "PhysicalForce.cuh"
#include "AllParticles.cuh"
#include <vector>
#include <memory>

class Simulation {
        public:

        Simulation(std::shared_ptr<AllParticles>& parts): parts(parts){};
        
        void simulate_one_step(const std::shared_ptr<AbstractForce>& force,const int num_parts,const double size);

        void init_particles(const double size, const int part_seed);

        void save_output(std::ofstream& fsave, int step);
        void save(std::ofstream& fsave);
        
        // private:
        std::shared_ptr<AllParticles> parts;


};

#endif