#include "common.cuh"
#include "AllParticles.cuh"

#include <memory>
#include <stdexcept>
#include <cmath>

#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <sm_60_atomic_functions.h>
#include <fstream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>
#include <thrust/execution_policy.h>






__global__ void move_kernel(double* dx, double* dy, double* dz,
                        double* dvx, double* dvy, double* dvz,
                        double* dax, double* day, double* daz, const double size, const int num_parts){
    // double size = dsize[0];
    // int num_parts = dnum_parts[0];
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i>=num_parts) return;
    dvx[i] += dax[i] * dt;
    dvy[i] += day[i] * dt;
    dvz[i] += daz[i] * dt;
    dx[i] += dvx[i] * dt;
    dy[i] += dvy[i] * dt;
    dz[i] += dvz[i] * dt;

    // Bounce from walls
    while (dx[i] < 0 || dx[i] > size) {
        dx[i] = (dx[i] < 0 ? -dx[i] : 2 * size - dx[i]);
        dvx[i] = -dvx[i];
    }

    while (dy[i] < 0 || dy[i] > size) {
        dy[i] = (dy[i] < 0 ? -dy[i] : 2 * size - dy[i]);
        dvy[i] = -dvy[i];
    }

    while (dz[i] < 0 || dz[i] > size) {
        dz[i] = (dz[i] < 0 ? -dz[i] : 2 * size - dz[i]);
        dvz[i] = -dvz[i];
    }
};



__global__ void ResetAcc(double* ax, double* ay, double* az, const int num_parts){
  unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=num_parts) return;
  ax[i] = 0.0;
  ay[i] = 0.0;
  az[i] = 0.0;
};

void AllParticles::ResetAccelerations(){
  ResetAcc<<<ceil(num_parts/th_per_block), th_per_block>>>(dax, day, daz, num_parts);
    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
      return;
    }
}

void AllParticles::init(){
    {
            thrust::default_random_engine rng;
            thrust::uniform_real_distribution<double> dist(0.0, size);
            thrust::uniform_real_distribution<double> dist1(-1.0, 1.0);

            for(int i=0; i<num_parts; i++){
              x_h[i] = dist(rng);
              y_h[i] = dist(rng);
              z_h[i] = dist(rng);
              vx[i] = dist1(rng);
              vy[i] = dist1(rng);
              vz[i] = dist1(rng);
              //pos[i] = make_double3(dist(rng), dist(rng), dist(rng));
              //vel[i] = make_double3(dist1(rng), dist1(rng), dist1(rng));
              //masses[i] = (dist1(rng) + 1.0);
              masses[i] = dist1(rng)+2.0;
              charges[i] = dist1(rng)*1e-19;
            }

            thrust::copy(x_h.begin(), x_h.end(), x.begin());
            thrust::copy(y_h.begin(), y_h.end(), y.begin());
            thrust::copy(z_h.begin(), z_h.end(), z.begin());

            // TODO mettere inizializzazione di xh e poi copy al vettore trust
            

            cudaDeviceSynchronize();
            ResetAccelerations();
            



        }
}


void AllParticles::save(std::ofstream& fsave){
    
  static bool first = true;

  if (first) {
      fsave << num_parts << " " << size << " " << nsteps << "\n";
      first = false;
  }
  //dovrei scrivere x_h[i] per risparmiare tempo ma non funziona. Ci penserò più tardi

  // Opzione 1:
  // thrust::copy(x.begin(), x.end(), x_h.begin());
  // thrust::copy(y.begin(), y.end(), y_h.begin());
  // thrust::copy(z.begin(), z.end(), z_h.begin());
  // Opzione 2:
  // x_h = x;
  // y_h = y;
  // z_h = z;
  cudaDeviceSynchronize();
  for(size_t i = 0; i < num_parts; i++){
    // TODO X_H
        fsave <<  x_h[i] << " " << y_h[i] << " " << z_h[i] << std::endl;
  }
  cudaDeviceSynchronize();
};

void AllParticles::save_output(std::ofstream& fsave, int step){
    // TODO FIX
    thrust::copy(x.begin(), x.end(), x_h.begin());
    thrust::copy(y.begin(), y.end(), y_h.begin());
    thrust::copy(z.begin(), z.end(), z_h.begin());
    save(fsave);
    cudaDeviceSynchronize();
    if(step > 0){
        if (step%10 == 0){
        fflush(stdout);
        printf("[ %d% ]\r", (int)(step*100/nsteps));
        }
    }
};

    // Integrate the ODE

    void AllParticles:: move() {
      move_kernel<<<grid_sizes, block_sizes>>> ( dx, dy, dz,
                                                 dvx, dvy, dvz,
                                                 dax, day, daz, size, num_parts);
                                                 
    };