#include "common.h"
#include "AllParticles.cuh"

#include <memory>
#include <stdexcept>
#include <cmath>


    // Integrate the ODE

    void AllParticles:: move(const double size, double* dx, double* dy, double* dz,
                                    double* dvx, double* dvy, double* dvz,
                                    double* dax, double* day, double* daz, const double size, const int num_parts , dim3 grid_sizes, const dim3 block_sizes) {
        move_kernel<<<grid_sizes, block_sizes>>> (double* dx, double* dy, double* dz,
                                                    double* dvx, double* dvy, double* dvz,
                                                    double* dax, double* day, double* daz, const double size, const int num_parts)
    };




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




