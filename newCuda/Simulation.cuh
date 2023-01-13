#ifndef HH__SIMULATION__HH
#define HH__SIMULATION__HH

#include "PhysicalForce.cuh"
#include "AllParticles.cuh"
#include <vector>
#include <memory>

class Simulation {
        public:

        Simulation(std::shared_ptr<AllParticles>& parts,  const int collision ): parts(parts), collision(collision){};
        
        void simulate_one_step(const std::shared_ptr<AbstractForce>& force,const int num_parts,const double size,  const int collision);

        void init_particles(const double size, const int part_seed);
        
        // private:
        std::shared_ptr<AllParticles> parts;
        int collision;


};

#endif