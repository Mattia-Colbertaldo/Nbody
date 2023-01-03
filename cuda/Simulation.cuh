#ifndef HH__SIMULATION__HH
#define HH__SIMULATION__HH

#include "PhysicalForce.cuh"
#include <vector>
#include <memory>
#include <cuda.h>

class Simulation {
        public:

        
        __host__ __device__ void simulate_one_step(AbstractForce* force,const int num_parts,const double size);

         __host__ __device__ void init_particles(const double size, const int part_seed);
        
        Particles D_parts;


};

#endif