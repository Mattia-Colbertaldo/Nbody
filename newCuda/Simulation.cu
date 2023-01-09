#include "AllParticles.cuh"
#include "Simulation.cuh"
#include <memory>
#include <stdexcept>
#include <cmath>
#include <random>
#include <iostream>


// Particle Initialization
void Simulation::init( const double size, const int part_seed) {
    /*
        input :  1. parts     : vettore di particelle
                 2. masses    : vettore delle masse
                 3. size      : dimensione della particella
                 4. part_seed : seme randomico
    */

    thrust::default_random_engine rng;
    thrust::uniform_real_distribution<double> dist(0.0, size);
    thrust::uniform_real_distribution<double> dist1(-1.0, 1.0);
    // TODO mettere inizializzazione di xh e poi copy al vettore trust
    for(int i=0; i<num_parts; i++){
      x[i] = dist(rng);
      y[i] = dist(rng);
      z[i] = dist(rng);
      // thrust::copy(x.begin(), x.end(), x_h.begin());
      // thrust::copy(y.begin(), y.end(), y_h.begin());
      // thrust::copy(z.begin(), z.end(), z_h.begin());
      vx[i] = dist1(rng);
      vy[i] = dist1(rng);
      vz[i] = dist1(rng);
      //pos[i] = make_double3(dist(rng), dist(rng), dist(rng));
      //vel[i] = make_double3(dist1(rng), dist1(rng), dist1(rng));
      //masses[i] = (dist1(rng) + 1.0);

      masses[i] = dist1(rng)+2.0;
      charges[i] = dist1(rng)*1e-19;
    }
    thrust::fill(acc.begin(), acc.end(), make_double3(0.0, 0.0, 0.0));
    parts.dx = thrust::raw_pointer_cast(x.data());
    parts.dy = thrust::raw_pointer_cast(y.data());
    parts.dz = thrust::raw_pointer_cast(z.data());
    parts.hx = thrust::raw_pointer_cast(x.data());
    parts.hy = thrust::raw_pointer_cast(y.data());
    parts.hz = thrust::raw_pointer_cast(z.data());
    parts.dvx = thrust::raw_pointer_cast(vx.data());
    parts.dvy = thrust::raw_pointer_cast(vy.data());
    parts.dvz = thrust::raw_pointer_cast(vz.data());
    parts.dax = thrust::raw_pointer_cast(ax.data());
    parts.day = thrust::raw_pointer_cast(ay.data());
    parts.daz = thrust::raw_pointer_cast(az.data());
    parts.dmasses = thrust::raw_pointer_cast(masses.data());
    parts.dcharges = thrust::raw_pointer_cast(charges.data());

    cudaDeviceSynchronize();
    ResetAcc<<<ceil(num_parts/th_per_block), th_per_block>>>(dax, day, daz, num_parts);
    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
      return;
    }



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
  th_per_block = fmin(32, num_parts);
  block_sizes.x = block_sizes.y = BLOCK_DIM;
  block_sizes.z = 1;
  grid_sizes.x = ceil(((double)num_parts)/((double)(block_sizes.x)));
  grid_sizes.y = ceil(((double)num_parts)/((double)(block_sizes.y)));
  grid_sizes.z = 1;
  if(first) std::cout << "GRID SIZE: " << grid_sizes.x << std::endl;
  ResetAcc<<<ceil((double)(num_parts)/(double)1024), 1024>>>(dax, day, daz, num_parts);
  force->force_application(dx, dy, dz, dvx, dvy, dvz, dax, day, daz, dmasses, dcharges, num_parts, grid_sizes, block_sizes);
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
  
  parts->move(dx, dy, dz, dvx, dvy, dvz, dax, day, daz, size, num_parts, grid_sizes, block_sizes);
  cudaDeviceSynchronize();
  error = cudaGetLastError();
  if (error != cudaSuccess) {
    fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    return;
  }
  
  
  if(first) std::cout << "Moving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  if(first) first = 0;
};
