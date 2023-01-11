#include "AllParticles.cuh"
#include "Simulation.cuh"
#include "common.cuh"
#include <memory>
#include <stdexcept>
#include <cmath>
#include <random>
#include <iostream>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <sm_60_atomic_functions.h>

#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/random.h>

#include <thrust/transform.h>
#include <thrust/execution_policy.h>

#include <thrust/zip_function.h>
#include <thrust/iterator/zip_iterator.h>
#include <fstream>


// Particle Initialization
void Simulation::init_particles( const double size, const int part_seed) {
    /*
        input :  1. parts     : vettore di particelle
                 2. masses    : vettore delle masse
                 3. size      : dimensione della particella
                 4. part_seed : seme randomico
    */

    parts->init();

    
};


void Simulation::save_output(std::ofstream& fsave, int step){
  parts->save_output(fsave, step);
}
void Simulation::save(std::ofstream& fsave){
  parts->save(fsave);
};

    



void Simulation::simulate_one_step( const std::shared_ptr<AbstractForce>& force,const int num_parts,const double size) 
{
  static bool first = 1;
  long t;
  if(first) t = clock();
  
  // long t = clock();
  // apply_force(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, sum_ax, sum_ay, sum_az);
  long t1;
  if(first) t1 = clock();
  parts->th_per_block = fmin(32, num_parts);
  parts->block_sizes.x = parts->block_sizes.y = BLOCK_DIM;
  parts->block_sizes.z = 1;
  parts->grid_sizes.x = ceil(((double)num_parts)/((double)(parts->block_sizes.x)));
  parts->grid_sizes.y = ceil(((double)num_parts)/((double)(parts->block_sizes.y)));
  parts->grid_sizes.z = 1;
  if(first) std::cout << "GRID SIZE: " << parts->grid_sizes.x << std::endl;
  parts->ResetAccelerations();
  force->force_application(parts->dx, parts->dy, parts->dz, parts->dvx, parts->dvy, parts->dvz, parts->dax, parts->day, 
                            parts->daz, parts->dmasses, parts->dcharges, num_parts, parts->grid_sizes, parts->block_sizes);
  // force_kernel<<<ceil((double)(num_parts)/(double)1024), 1024>>>(dx, dy, dz, dax, day, daz, dmasses, dcharges, num_parts);
  if(first) std::cout << "Applying force: Kernel Loop: " << ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  
  cudaDeviceSynchronize();
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    return;
  };
  
  if(first) std::cout << "Applying force: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  if(first) t = clock();
  // <<<grid_dim, block_size>>>
  
  parts->move();
  cudaDeviceSynchronize();
  error = cudaGetLastError();
  if (error != cudaSuccess) {
    fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    return;
  }
  
  
  if(first) std::cout << "Moving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  if(first) first = 0;
};


