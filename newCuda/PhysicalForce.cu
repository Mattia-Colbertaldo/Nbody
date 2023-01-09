#include "common.cuh"
#include "PhysicalForce.cuh"
#include <cuda.h>

    
__global__ void RepulsiveForce :: 
kernel_no_tiling_force(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts){
    int thx = threadIdx.x + blockDim.x * blockIdx.x;
    int thy = threadIdx.y + blockDim.y * blockIdx.y;
    

    // printf("%d, %d\n", thx, thy);
    // se io sono il thread (3,4) applico la forza a 3 e a 4
    // lo faccio solo per i thread la cui x < y
    if(thx < thy && thy < num_parts){
      double dx = x[thy] - x[thx];
      double dy = y[thy] - y[thx];
      double dz = z[thy] - z[thx];
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > cutoff * cutoff) return;
      // *** EXPERIMENTAL *** //
      if(r2 < min_r*min_r){
        
        // spingo l'altra particella : applico alla mia vicina una forza uguale a F = m * a ,
        // quindi applico a lei un'accelerazione di - a_mia * m_mia/m_sua
        // atomicAdd((double*)(ax + thy), (double)ax[thx]*masses[thx]/masses[thy]);
        // atomicAdd((double*)(ax + thx), (double)-ax[thy]*masses[thy]/masses[thx]);

        // atomicAdd((double*)(ay + thy), (double)ax[thx]*masses[thx]/masses[thy]);
        // atomicAdd((double*)(ay + thx), (double)-ax[thy]*masses[thy]/masses[thx]);

        // atomicAdd((double*)(az + thy), (double)ax[thx]*masses[thx]/masses[thy]);
        // atomicAdd((double*)(az + thx), (double)-ax[thy]*masses[thy]/masses[thx]);

        double mx = masses[thx];
        double my = masses[thy];
        // TODO ARGUMENT
        if(1){
          // URTO ANELASTICO:
          vx[thx] = (double)(mx*vx[thx] + my*vx[thy])/(double)(mx+my);
          vx[thy] = (double)(my*vx[thy] + mx*vx[thx])/(double)(my+mx);

          vy[thx] = (double)(mx*vy[thx] + my*vy[thy])/(double)(mx+my);
          vy[thy] = (double)(my*vy[thy] + mx*vy[thx])/(double)(my+mx);

          vz[thx] = (double)(mx*vz[thx] + my*vz[thy])/(double)(mx+my);
          vz[thy] = (double)(my*vz[thy] + mx*vz[thx])/(double)(my+mx);
        }
        else{
          // URTO ELASTICO
          vx[thx] = (double)vx[thx]*(mx-my)/(mx + my) + 2*vx[thy]*my/(mx+my);
          vx[thy] = (double)vx[thy]*(my-mx)/(mx + my) + 2*vx[thx]*mx/(mx+my);

          vy[thx] = (double)vy[thx]*(mx-my)/(mx + my) + 2*vy[thy]*my/(mx+my);
          vy[thy] = (double)vy[thy]*(my-mx)/(mx + my) + 2*vy[thx]*mx/(mx+my);

          vz[thx] = (double)vz[thx]*(mx-my)/(mx + my) + 2*vz[thy]*my/(mx+my);
          vz[thy] = (double)vz[thy]*(my-mx)/(mx + my) + 2*vz[thx]*mx/(mx+my);
        }
        return;
      }
      // ******************** //
      r2 = fmax(r2, min_r * min_r);
      double coef = G / r2;

      atomicAdd((double*)(ax + thx), (double)coef*dx*masses[thy]);
      atomicAdd((double*)(ax + thy), (double)-coef*dx*masses[thx]);
      // ax[index] += coef*dx;
      atomicAdd((double*)(ay + thx), (double)coef*dy*masses[thy]);
      atomicAdd((double*)(ay + thy), (double)-coef*dy*masses[thx]);
      // ay[index] += coef*dy;
      atomicAdd((double*)(az + thx), (double)coef*dz*masses[thy]);
      atomicAdd((double*)(az + thy), (double)-coef*dz*masses[thx]);
      // az[index] += coef*dz;

    }
    __syncthreads();
  
}

 /* 
void GravitationalForce :: force_application() const {
    Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef =  (G * neighbor.mass / r2) ;

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;

    // neighbor.ax += coef * dx;
    // neighbor.ay += coef * dy;
    // neighbor.az += coef * dz;
};


void GravitationalAssistForce:: force_application() const {
    Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef;

    // Very simple short-range repulsive force
    if(r2>0.0001){
        coef =  G * neighbor.mass / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        coef = -( G * neighbor.mass / r2 ) * 3 ;
    }

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
};



     
void ProtonForce :: force_application() const {
    Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef =  K * proton_charge * proton_charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
};


   
void CoulombForce :: force_application() const {
    
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef = std::pow(scale, 2) * K * p.charge * neighbor.charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
};




*/