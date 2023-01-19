#include "common.cuh"
#include "PhysicalForce.cuh"
#include <cuda.h>
  

using namespace common_h;

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

__global__ void
kernel_super_tiling_force_gravitational(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
                        const int collision ){

    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;    

    const int tile_size = common_block_size;
    __shared__ double tile1x[tile_size];
    __shared__ double tile1y[tile_size];
    __shared__ double tile1z[tile_size];
    __shared__ double tile1_masses[tile_size];
    __shared__ double tile1_charges[tile_size];
    __shared__ double tile1_vx[tile_size];
    __shared__ double tile1_vy[tile_size];
    __shared__ double tile1_vz[tile_size];

    __shared__ double tile2x[tile_size];
    __shared__ double tile2y[tile_size];
    __shared__ double tile2z[tile_size];
    __shared__ double tile2_masses[tile_size];
    __shared__ double tile2_charges[tile_size];
    __shared__ double tile2_vx[tile_size];
    __shared__ double tile2_vy[tile_size];
    __shared__ double tile2_vz[tile_size];

    for(int i=0; i<tile_size; i++){
      if(threadIdx.x < tile_size && threadIdx.y < tile_size && thy < num_parts && thx < num_parts && thx != thy ){
        tile1x[threadIdx.x] = x[thx];
        tile1y[threadIdx.x] = y[thx];
        tile1z[threadIdx.x] = z[thx];
        tile1_masses[threadIdx.x] = masses[thx];
        tile1_vx[threadIdx.x] = vx[threadIdx.x];
        tile1_vy[threadIdx.x] = vy[threadIdx.x];
        tile1_vz[threadIdx.x] = vz[threadIdx.x];

        tile2x[threadIdx.y] = x[thy];
        tile2y[threadIdx.y] = y[thy];
        tile2z[threadIdx.y] = z[thy];
        tile2_masses[threadIdx.y] = masses[thy];
        tile2_vx[threadIdx.y] = vx[threadIdx.y];
        tile2_vy[threadIdx.y] = vy[threadIdx.y];
        tile2_vz[threadIdx.y] = vz[threadIdx.y];
      }
    }
    

    
    double mx = tile1_masses[threadIdx.x];
    double my = tile2_masses[threadIdx.y];

    if(thx < num_parts && thy < num_parts && thx!=thy){
      double dx = tile2x[threadIdx.y] - tile1x[threadIdx.x];
      double dy = tile2y[threadIdx.y] - tile1y[threadIdx.x];
      double dz = tile2z[threadIdx.y] - tile1z[threadIdx.x];
      double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
      if (r2 > cutoff * cutoff) return;
      
      if(r2 < min_r*min_r){
        if(collision > 0){
            
        if(collision== 1){
          // inelastic collision
          vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);

          vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);

          vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
        }
        
        else if(collision== 2){
          // elastic collision
          vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);

          vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);

          vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
        }
        }
        return;
      }
      
      r2 = fmax(r2, min_r * min_r);
      double coef =  (G / r2) ;

      atomicAdd((double*)(ax + thx), (double)coef*dx*my);
      
      atomicAdd((double*)(ay + thx), (double)coef*dy*my);
      
      atomicAdd((double*)(az + thx), (double)coef*dz*my);
      

    }
  
}

__global__ void
kernel_tiling_force_repulsive(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
                        const int collision ){
    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;
    if(thx >= thy) return;
    const int tile_size = common_block_size;
    __shared__ double tile1x[tile_size];
    __shared__ double tile1y[tile_size];
    __shared__ double tile1z[tile_size];
    __shared__ double tile1_masses[tile_size];
    __shared__ double tile1_charges[tile_size];
    __shared__ double tile1_vx[tile_size];
    __shared__ double tile1_vy[tile_size];
    __shared__ double tile1_vz[tile_size];

    __shared__ double tile2x[tile_size];
    __shared__ double tile2y[tile_size];
    __shared__ double tile2z[tile_size];
    __shared__ double tile2_masses[tile_size];
    __shared__ double tile2_charges[tile_size];
    __shared__ double tile2_vx[tile_size];
    __shared__ double tile2_vy[tile_size];
    __shared__ double tile2_vz[tile_size];

    for(int i=0; i<tile_size; i++){
      if(threadIdx.x < tile_size && threadIdx.y < tile_size && thy < num_parts && thx < thy ){
        tile1x[threadIdx.x] = x[thx];
        tile1y[threadIdx.x] = y[thx];
        tile1z[threadIdx.x] = z[thx];
        tile1_masses[threadIdx.x] = masses[thx];
        tile1_vx[threadIdx.x] = vx[threadIdx.x];
        tile1_vy[threadIdx.x] = vy[threadIdx.x];
        tile1_vz[threadIdx.x] = vz[threadIdx.x];

        tile2x[threadIdx.y] = x[thy];
        tile2y[threadIdx.y] = y[thy];
        tile2z[threadIdx.y] = z[thy];
        tile2_masses[threadIdx.y] = masses[thy];
        tile2_vx[threadIdx.y] = vx[threadIdx.y];
        tile2_vy[threadIdx.y] = vy[threadIdx.y];
        tile2_vz[threadIdx.y] = vz[threadIdx.y];
      }
    }
    

    
    double mx = tile1_masses[threadIdx.x];
    double my = tile2_masses[threadIdx.y];

    
    if(thx < thy && thy < num_parts){
      double dx = tile2x[threadIdx.y] - tile1x[threadIdx.x];
      double dy = tile2y[threadIdx.y] - tile1y[threadIdx.x];
      double dz = tile2z[threadIdx.y] - tile1z[threadIdx.x];
      double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
      if (r2 > cutoff * cutoff) return;
      
      if(r2 < min_r*min_r){
        if(collision > 0){
        
        if(collision== 1){
          // inelastic collision
          vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
          vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

          vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
          vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

          vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
          vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
        }
        
        else if(collision== 2){
          // elastic collision
          vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
          vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

          vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
          vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

          vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
          vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
        }
        }
        return;
      }
      r2 = fmax(r2, min_r * min_r);
      double r = std::sqrt(r2);
      double coef = 1e-2*(1 - cutoff / r) / r2 ;

      atomicAdd((double*)(ax + thx), (double)coef*dx*masses[thy]);
      atomicAdd((double*)(ax + thy), (double)-coef*dx*masses[thx]);
      
      atomicAdd((double*)(ay + thx), (double)coef*dy*masses[thy]);
      atomicAdd((double*)(ay + thy), (double)-coef*dy*masses[thx]);
     
      atomicAdd((double*)(az + thx), (double)coef*dz*masses[thy]);
      atomicAdd((double*)(az + thy), (double)-coef*dz*masses[thx]);
      

    }
    
  
}

__global__ void
kernel_tiling_force_gravitational(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
                        const int collision ){

    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;
    if(thx >= thy) return;
    

    const int tile_size = common_block_size;
    __shared__ double tile1x[tile_size];
    __shared__ double tile1y[tile_size];
    __shared__ double tile1z[tile_size];
    __shared__ double tile1_masses[tile_size];
    __shared__ double tile1_charges[tile_size];
    __shared__ double tile1_vx[tile_size];
    __shared__ double tile1_vy[tile_size];
    __shared__ double tile1_vz[tile_size];

    __shared__ double tile2x[tile_size];
    __shared__ double tile2y[tile_size];
    __shared__ double tile2z[tile_size];
    __shared__ double tile2_masses[tile_size];
    __shared__ double tile2_charges[tile_size];
    __shared__ double tile2_vx[tile_size];
    __shared__ double tile2_vy[tile_size];
    __shared__ double tile2_vz[tile_size];

    for(int i=0; i<tile_size; i++){
      if(threadIdx.x < tile_size && threadIdx.y < tile_size && thy < num_parts && thx < thy ){
        tile1x[threadIdx.x] = x[thx];
        tile1y[threadIdx.x] = y[thx];
        tile1z[threadIdx.x] = z[thx];
        tile1_masses[threadIdx.x] = masses[thx];
        tile1_vx[threadIdx.x] = vx[threadIdx.x];
        tile1_vy[threadIdx.x] = vy[threadIdx.x];
        tile1_vz[threadIdx.x] = vz[threadIdx.x];

        tile2x[threadIdx.y] = x[thy];
        tile2y[threadIdx.y] = y[thy];
        tile2z[threadIdx.y] = z[thy];
        tile2_masses[threadIdx.y] = masses[thy];
        tile2_vx[threadIdx.y] = vx[threadIdx.y];
        tile2_vy[threadIdx.y] = vy[threadIdx.y];
        tile2_vz[threadIdx.y] = vz[threadIdx.y];
      }
    }
    

    
    double mx = tile1_masses[threadIdx.x];
    double my = tile2_masses[threadIdx.y];

    
    if(thx < thy && thy < num_parts){
      double dx = tile2x[threadIdx.y] - tile1x[threadIdx.x];
      double dy = tile2y[threadIdx.y] - tile1y[threadIdx.x];
      double dz = tile2z[threadIdx.y] - tile1z[threadIdx.x];
      double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
      if (r2 > cutoff * cutoff) return;
      
      if(r2 < min_r*min_r){
        if(collision > 0){
      
        if(collision== 1){
          // inelastic collision
          vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
          vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

          vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
          vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

          vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
          vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
        }
        
        else if(collision== 2){
          // elastic collision
          vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
          vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

          vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
          vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

          vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
          vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
        }
        }
        return;
      }
      r2 = fmax(r2, min_r * min_r);
      double coef =  (G / r2) ;

      atomicAdd((double*)(ax + thx), (double)coef*dx*my);
      atomicAdd((double*)(ax + thy), (double)-coef*dx*mx);
    
      atomicAdd((double*)(ay + thx), (double)coef*dy*my);
      atomicAdd((double*)(ay + thy), (double)-coef*dy*mx);
    
      atomicAdd((double*)(az + thx), (double)coef*dz*my);
      atomicAdd((double*)(az + thy), (double)-coef*dz*mx);
      

    }
    
  
}



__global__ void
kernel_tiling_force_gravitational_assist(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
                        const int collision ){
    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;
    if(thx >= thy) return;
    const int tile_size = common_block_size;
    __shared__ double tile1x[tile_size];
    __shared__ double tile1y[tile_size];
    __shared__ double tile1z[tile_size];
    __shared__ double tile1_masses[tile_size];
    __shared__ double tile1_charges[tile_size];
    __shared__ double tile1_vx[tile_size];
    __shared__ double tile1_vy[tile_size];
    __shared__ double tile1_vz[tile_size];

    __shared__ double tile2x[tile_size];
    __shared__ double tile2y[tile_size];
    __shared__ double tile2z[tile_size];
    __shared__ double tile2_masses[tile_size];
    __shared__ double tile2_charges[tile_size];
    __shared__ double tile2_vx[tile_size];
    __shared__ double tile2_vy[tile_size];
    __shared__ double tile2_vz[tile_size];

    for(int i=0; i<tile_size; i++){
      if(threadIdx.x < tile_size && threadIdx.y < tile_size && thy < num_parts && thx < thy ){
        tile1x[threadIdx.x] = x[thx];
        tile1y[threadIdx.x] = y[thx];
        tile1z[threadIdx.x] = z[thx];
        tile1_masses[threadIdx.x] = masses[thx];
        tile1_vx[threadIdx.x] = vx[threadIdx.x];
        tile1_vy[threadIdx.x] = vy[threadIdx.x];
        tile1_vz[threadIdx.x] = vz[threadIdx.x];

        tile2x[threadIdx.y] = x[thy];
        tile2y[threadIdx.y] = y[thy];
        tile2z[threadIdx.y] = z[thy];
        tile2_masses[threadIdx.y] = masses[thy];
        tile2_vx[threadIdx.y] = vx[threadIdx.y];
        tile2_vy[threadIdx.y] = vy[threadIdx.y];
        tile2_vz[threadIdx.y] = vz[threadIdx.y];
      }
    }
    

    
    double mx = tile1_masses[threadIdx.x];
    double my = tile2_masses[threadIdx.y];

    
    if(thx < thy && thy < num_parts){
      double dx = tile2x[threadIdx.y] - tile1x[threadIdx.x];
      double dy = tile2y[threadIdx.y] - tile1y[threadIdx.x];
      double dz = tile2z[threadIdx.y] - tile1z[threadIdx.x];
      double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
      if (r2 > cutoff * cutoff) return;
      
      if(r2 < min_r*min_r){
        if(collision > 0){
        
        if(collision== 1){
          // inelastic collision
          vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
          vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

          vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
          vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

          vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
          vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
        }
        
        else if(collision== 2){
          // elastic collision
          vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
          vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

          vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
          vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

          vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
          vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
        }
        }
        return;
      }
      r2 = fmax(r2, min_r * min_r);
      double coef;

        // Very simple short-range repulsive force
        if(r2>0.0001){
            coef =  G  / r2 ;
        }
        else
        //gravity-assist : repulsive force
        {
            coef = -( G  / r2 ) * 3 ;
        }

      atomicAdd((double*)(ax + thx), (double)coef*dx*masses[thy]);
      atomicAdd((double*)(ax + thy), (double)-coef*dx*masses[thx]);
      
      atomicAdd((double*)(ay + thx), (double)coef*dy*masses[thy]);
      atomicAdd((double*)(ay + thy), (double)-coef*dy*masses[thx]);
      
      atomicAdd((double*)(az + thx), (double)coef*dz*masses[thy]);
      atomicAdd((double*)(az + thy), (double)-coef*dz*masses[thx]);
     

    }
    
  
}


__global__ void
kernel_tiling_force_proton(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
                        const int collision ){
    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;
    if(thx >= thy) return;
    const int tile_size = common_block_size;
    __shared__ double tile1x[tile_size];
    __shared__ double tile1y[tile_size];
    __shared__ double tile1z[tile_size];
    __shared__ double tile1_masses[tile_size];
    __shared__ double tile1_charges[tile_size];
    __shared__ double tile1_vx[tile_size];
    __shared__ double tile1_vy[tile_size];
    __shared__ double tile1_vz[tile_size];

    __shared__ double tile2x[tile_size];
    __shared__ double tile2y[tile_size];
    __shared__ double tile2z[tile_size];
    __shared__ double tile2_masses[tile_size];
    __shared__ double tile2_charges[tile_size];
    __shared__ double tile2_vx[tile_size];
    __shared__ double tile2_vy[tile_size];
    __shared__ double tile2_vz[tile_size];

    for(int i=0; i<tile_size; i++){
      if(threadIdx.x < tile_size && threadIdx.y < tile_size && thy < num_parts && thx < thy ){
        tile1x[threadIdx.x] = x[thx];
        tile1y[threadIdx.x] = y[thx];
        tile1z[threadIdx.x] = z[thx];
        tile1_masses[threadIdx.x] = masses[thx];
        tile1_vx[threadIdx.x] = vx[threadIdx.x];
        tile1_vy[threadIdx.x] = vy[threadIdx.x];
        tile1_vz[threadIdx.x] = vz[threadIdx.x];

        tile2x[threadIdx.y] = x[thy];
        tile2y[threadIdx.y] = y[thy];
        tile2z[threadIdx.y] = z[thy];
        tile2_masses[threadIdx.y] = masses[thy];
        tile2_vx[threadIdx.y] = vx[threadIdx.y];
        tile2_vy[threadIdx.y] = vy[threadIdx.y];
        tile2_vz[threadIdx.y] = vz[threadIdx.y];
      }
    }
    

    
    double mx = tile1_masses[threadIdx.x];
    double my = tile2_masses[threadIdx.y];

    
    if(thx < thy && thy < num_parts){
      double dx = tile2x[threadIdx.y] - tile1x[threadIdx.x];
      double dy = tile2y[threadIdx.y] - tile1y[threadIdx.x];
      double dz = tile2z[threadIdx.y] - tile1z[threadIdx.x];
      double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
      if (r2 > cutoff * cutoff) return;
      
      if(r2 < min_r*min_r){
        if(collision > 0){
        
        
        if(collision== 1){
          // inelastic collision
          vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
          vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

          vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
          vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

          vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
          vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
        }
        
        else if(collision== 2){
          // elastic collision
          vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
          vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

          vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
          vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

          vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
          vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
          }
        }
        return;
      }
      // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef =  K * proton_charge * proton_charge / r2  ;
    
      
      double charge_product = charges[thy]*charges[thx] ;
      atomicAdd((double*)(ax + thx), (double)coef*dx*+charge_product);
      atomicAdd((double*)(ax + thy), (double)-coef*dx*+charge_product);
      
      atomicAdd((double*)(ay + thx), (double)coef*dy*+charge_product);
      atomicAdd((double*)(ay + thy), (double)-coef*dy*+charge_product);
     
      atomicAdd((double*)(az + thx), (double)coef*dz*+charge_product);
      atomicAdd((double*)(az + thy), (double)-coef*dz*+charge_product);
      

    }
    
  
}



__global__ void
kernel_tiling_force_coulomb(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
                        const int collision ){
    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;
    if(thx >= thy) return;
    const int tile_size = common_block_size;
    __shared__ double tile1x[tile_size];
    __shared__ double tile1y[tile_size];
    __shared__ double tile1z[tile_size];
    __shared__ double tile1_masses[tile_size];
    __shared__ double tile1_charges[tile_size];
    __shared__ double tile1_vx[tile_size];
    __shared__ double tile1_vy[tile_size];
    __shared__ double tile1_vz[tile_size];

    __shared__ double tile2x[tile_size];
    __shared__ double tile2y[tile_size];
    __shared__ double tile2z[tile_size];
    __shared__ double tile2_masses[tile_size];
    __shared__ double tile2_charges[tile_size];
    __shared__ double tile2_vx[tile_size];
    __shared__ double tile2_vy[tile_size];
    __shared__ double tile2_vz[tile_size];

    for(int i=0; i<tile_size; i++){
      if(threadIdx.x < tile_size && threadIdx.y < tile_size && thy < num_parts && thx < thy ){
        tile1x[threadIdx.x] = x[thx];
        tile1y[threadIdx.x] = y[thx];
        tile1z[threadIdx.x] = z[thx];
        tile1_masses[threadIdx.x] = masses[thx];
        tile1_vx[threadIdx.x] = vx[threadIdx.x];
        tile1_vy[threadIdx.x] = vy[threadIdx.x];
        tile1_vz[threadIdx.x] = vz[threadIdx.x];

        tile2x[threadIdx.y] = x[thy];
        tile2y[threadIdx.y] = y[thy];
        tile2z[threadIdx.y] = z[thy];
        tile2_masses[threadIdx.y] = masses[thy];
        tile2_vx[threadIdx.y] = vx[threadIdx.y];
        tile2_vy[threadIdx.y] = vy[threadIdx.y];
        tile2_vz[threadIdx.y] = vz[threadIdx.y];
      }
    }
    

    
    double mx = tile1_masses[threadIdx.x];
    double my = tile2_masses[threadIdx.y];

  

    
    if(thx < thy && thy < num_parts){
      double dx = tile2x[threadIdx.y] - tile1x[threadIdx.x];
      double dy = tile2y[threadIdx.y] - tile1y[threadIdx.x];
      double dz = tile2z[threadIdx.y] - tile1z[threadIdx.x];
      double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
      if (r2 > cutoff * cutoff) return;
      
      if(r2 < min_r*min_r){
        if(collision > 0){
        
        
        if(collision== 1){
          // inelastic collision
          vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
          vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

          vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
          vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

          vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
          vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
        }
        
        else if(collision== 2){
          // elastic collision
          vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
          vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

          vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
          vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

          vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
          vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
        }
        }
        return;
      }
      
      r2 = fmax(r2, min_r * min_r);
      double coef = std::pow(scale, 2) * K  / r2  ;
      double charge_product = charges[thy]*charges[thx] ;
      atomicAdd((double*)(ax + thx), (double)coef*dx*+charge_product);
      atomicAdd((double*)(ax + thy), (double)-coef*dx*+charge_product);
    
      atomicAdd((double*)(ay + thx), (double)coef*dy*+charge_product);
      atomicAdd((double*)(ay + thy), (double)-coef*dy*+charge_product);
      
      atomicAdd((double*)(az + thx), (double)coef*dz*+charge_product);
      atomicAdd((double*)(az + thy), (double)-coef*dz*+charge_product);
    

    }
    
  
}

void RepulsiveForce :: force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
    double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
    const int collision , dim3 grid_sizes, const dim3 block_sizes ) const {
    
    kernel_tiling_force_repulsive<<<grid_sizes, block_sizes>>>(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, collision);
}

void GravitationalForce :: force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
    double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
    const int collision , dim3 grid_sizes, const dim3 block_sizes ) const {
    
    kernel_super_tiling_force_gravitational<<<grid_sizes, block_sizes>>>(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, collision);
}


void GravitationalAssistForce :: force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
    double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
    const int collision , dim3 grid_sizes, const dim3 block_sizes ) const {
    
    kernel_tiling_force_gravitational_assist<<<grid_sizes, block_sizes>>>(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, collision);
}


void ProtonForce :: force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
    double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
    const int collision , dim3 grid_sizes, const dim3 block_sizes ) const {
    
    kernel_tiling_force_proton<<<grid_sizes, block_sizes>>>(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, collision);
}

void CoulombForce :: force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
    double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, 
    const int collision , dim3 grid_sizes, const dim3 block_sizes ) const {
    
    kernel_tiling_force_coulomb<<<grid_sizes, block_sizes>>>(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, collision);
}
