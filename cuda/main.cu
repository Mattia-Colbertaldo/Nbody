// Tenere Force.cu
// Vedere se funziona all_particles e Simulatioin
// Aggiungere Arg e Output (che adesso sono in particles)

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
#include <fstream>

#include <thrust/transform.h>
#include <thrust/execution_policy.h>

#include "Output.cuh"
#include "Find_Arg.cuh"
#include "common.cuh"
#include "PhysicalForce.cuh"
#include "AllParticles.cuh"
#include "Simulation.cuh"

#include <time.h>


#include <thrust/zip_function.h>
#include <thrust/iterator/zip_iterator.h>



#include <cassert>
#include <iostream>






void DisplayHeader()
{
    const int kb = 1024;
    const int mb = kb * kb;
    std::cout << "NBody.GPU" << std::endl << "=========" << std::endl << std::endl;

    std::cout << "CUDA version:   v" << CUDART_VERSION << std::endl;

    int devCount;
    cudaGetDeviceCount(&devCount);
    std::cout << "CUDA Devices: " << std::endl << std::endl;

    for(int i = 0; i < devCount; ++i)
    {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        std::cout << i << ": " << props.name << ": " << props.major << "." << props.minor << std::endl;
        std::cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << std::endl;
        std::cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << std::endl;
        std::cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << std::endl;
        std::cout << "  Block registers: " << props.regsPerBlock << std::endl << std::endl;

        std::cout << "  Warp size:         " << props.warpSize << std::endl;
        std::cout << "  Threads per block: " << props.maxThreadsPerBlock << std::endl;
        std::cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << std::endl;
        std::cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << std::endl;
        std::cout << std::endl;
    }
}









/*

    // TODO: SPOSTARE LE FORZE DENTRO PHYSICALFORCE.CU E POI CHIAMARE FORCE.APPLY_FORCE()
// (NON E' DETTO CHE FUNZIONI, DOBBIAMO VEDERE SE LE CHIAMATE A FUNZIONI GLOBAL VANNO, FACCIAMOLO INSIEME CHE CI AVEVO GIA' PROVATO)
__global__ void kernel_test_force(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts){
    
    int thx = threadIdx.x;
    int thy = threadIdx.y;
    int Row = threadIdx.x + blockDim.x * blockIdx.x;
    int Col = threadIdx.y + blockDim.y * blockIdx.y;


    __shared__ double tilex_1[TILE_WIDTH];
    __shared__ double tiley_1[TILE_WIDTH];
    __shared__ double tilez_1[TILE_WIDTH];
    __shared__ double tile_masses_1[TILE_WIDTH];
    __shared__ double tile_charges_1[TILE_WIDTH];
    

    __shared__ double tilex_2[TILE_WIDTH];
    __shared__ double tiley_2[TILE_WIDTH];
    __shared__ double tilez_2[TILE_WIDTH];
    __shared__ double tile_masses_2[TILE_WIDTH];
    __shared__ double tile_charges_2[TILE_WIDTH];;
    
    // I do calculations only if x index is a successor of y index

    // Loop over all particles with tiles
    if(Row < 0 || Row >= num_parts) return;
    if(Col < 0 || Col >= num_parts) return;
    if(thx>TILE_WIDTH || thy>TILE_WIDTH) return;
    
    tilex_1[thx] = x[ Row ];
    tiley_1[thx] = y[ Row ];
    tilez_1[thx] = z[ Row ];
    tile_masses_1[thx] = masses[ Row ];
    tile_charges_1[thx] = charges[ Row ];
  
    tilex_2[thy] = x[ Col ];
    tiley_2[thy] = y[ Col ];
    tilez_2[thy] = z[ Col ];
    tile_masses_2[thy] = masses[ Col ];
    tile_charges_2[thy] = charges[ Col ];

    __syncthreads();

    if(Row<Col && thx<TILE_WIDTH && thy<TILE_WIDTH && Col < num_parts){
      double dx = tilex_2[ thy ] - tilex_1[ thx ];
      double dy = tiley_2[ thy ] - tiley_1[ thx ];
      double dz = tilez_2[ thy ] - tilez_1[ thx ];
      double r2 = dx * dx + dy * dy + dz * dz;

      
      double mx = tile_masses_1[thx];
      double my = tile_masses_2[thy];

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

        // TODO: CHANGE LOADING X[] INTO TILEX[]
        // TODO ARGUMENT
        if(0){
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

      atomicAdd((double*)(ax + Row), (double)coef*dx*my);
      atomicAdd((double*)(ax + Col), (double)-coef*dx*mx);
      // ax[index] += coef*dx;
      atomicAdd((double*)(ay + Row), (double)coef*dy*my);
      atomicAdd((double*)(ay + Col), (double)-coef*dy*mx);
      // ay[index] += coef*dy;
      atomicAdd((double*)(az + Row), (double)coef*dz*my);
      atomicAdd((double*)(az + Col), (double)-coef*dz*mx);
      // az[index] += coef*dz;

    }
    __syncthreads();
}

*/

  





int main(int argc, char** argv)
{

  DisplayHeader();

  Find_Arg finder= Find_Arg(argc, argv);
  if (finder.find_int_arg("-h", 0) >= 0) {
      std::cout << "Options:" << std::endl;
      std::cout << "-h: see this help" << std::endl;
      std::cout << "-n <int>: set number of particles" << std::endl;
      std::cout << "-o <filename>: set the output file name" << std::endl;
      std::cout << "-s <int>: set particle initialization seed" << std::endl;
      std::cout << "-t <int>: set number of threads (working only in parallel mode) [default = 8]" << std::endl;
      std::cout << "-f <int>: set force: default, repulsive, gravitational, assist, proton, coulomb" << std::endl;
      return 0;
  }

  // Open Output File
  std::string savename = finder.find_string_arg("-o", "out.txt");
  if (savename != "") std::cout << "Creating file " << savename << "..." << std::endl;
  
  //Find force
  std::string forcename = finder.find_string_arg("-f", "repulsive");
  if (forcename != "") std::cout << "Choosing " <<  forcename << " force..." << std::endl;
  else{
      std::string def="default";
      forcename= &def[0];
      std::cout << "Choosing default force..." << std::endl;;
  }

  std::shared_ptr<AbstractForce> force =finder.find_force(forcename) ;

  //find collision type
  int collision = finder.find_int_arg("-c", 0);
  std::cout << "Choosing " <<  collision << " collision type..." << std::endl;


  const int num_parts = finder.find_int_arg("-n", 1000);
  std::cout << "Starting simulation with " << num_parts << " particles." <<std::endl;
  const int part_seed = finder.find_int_arg("-s", 0);
  
  std::unique_ptr<AllParticles> all_particles = std::make_unique<AllParticles>(num_parts, force);
  
  const double size = std::sqrt(density * num_parts);
  const int num_th = finder.find_int_arg("-t", 8);
  long t = clock();
  Simulation s = Simulation(all_particles, collision );
  std::cout << "Initialization: ";
  s.init_particles(size, part_seed);
  cudaDeviceSynchronize();
  

  std::cout << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  // std::cout << "Moving data from GPU to CPU: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  t = clock();
  // thrust::copy(x.begin(), x.end(), x_h.begin());
  // thrust::copy(y.begin(), y.end(), y_h.begin());
  // thrust::copy(z.begin(), z.end(), z_h.begin());
  Output output= Output(num_parts, savename);
  output.save( s.parts, size, nsteps);
  std::cout << "Saving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  std::cout << "Now entering the for loop." << std::endl;
  t = clock();
  long t1;
  long savetime = 0.;
  for(int step=0; step<nsteps; step++){
    if(step == nsteps/2) t1 = clock();
    s.simulate_one_step(force, num_parts, size);
    cudaDeviceSynchronize();
    if(step == nsteps/2) std::cout << "Simulating one step: " << ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
    if(step == nsteps/2) t1 = clock();
    output.save_output( savefreq, s.parts , step, nsteps, size);
    if(step == nsteps/2) std::cout << "Saving to device buffer: " << ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
    else if(step == nsteps - 1){
      savetime = ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC;
      std::cout << "Final Saving: " << savetime << std::endl;
    }
    if(step == nsteps/2) std::cout << "One loop iteration: " << savetime << " ms" << std::endl;
  }
  long all_steps = ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC;
  std::cout << std::endl << "***************************************" << std::endl;
  std::cout << " + Calculations: " << all_steps-savetime << " ms [ " << ((all_steps-savetime)*100/all_steps) << " % ]" << std::endl;
  std::cout << " + Saving: " << savetime << " ms [ " << (savetime*100/all_steps) << " % ]" << std::endl;
  std::cout << " = Simulating all the steps: " << all_steps << " ms" << std::endl;
  std::cout << "***************************************" << std::endl << std::endl;;
  


  return 0;
}
