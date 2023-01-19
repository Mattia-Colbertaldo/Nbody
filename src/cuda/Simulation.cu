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


// Particle Initialization, 
// we call the method of AllParticles class so that it can initialize elements that are part of the class itself
void Simulation::init_particles( const double size, const int part_seed) {
    /*
        input :  1. parts     : vettore di particelle
                 2. masses    : vettore delle masse
                 3. size      : dimensione della particella
                 4. part_seed : seme randomico
    */

    parts->init();

    
};

    



void Simulation::simulate_one_step( const std::shared_ptr<AbstractForce>& force,const int num_parts,const double size) 
{
  static bool first = 1;
  long t;
  if(first) t = clock();
  long t1;
  if(first) t1 = clock();  
  if(first) std::cout << "GRID SIZE: " << parts->grid_sizes.x << " x " << parts->grid_sizes.y << std::endl;
  parts->ResetAccelerations();
  
  force->force_application(parts->dx, parts->dy, parts->dz, parts->dvx, parts->dvy, parts->dvz, parts->dax, parts->day, 
                            parts->daz, parts->dmasses, parts->dcharges, num_parts, collision, parts->grid_sizes, parts->block_sizes);
  
  cudaDeviceSynchronize();
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    return;
  };
  
  if(first) std::cout << "Applying force: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  if(first) t = clock();
  
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


