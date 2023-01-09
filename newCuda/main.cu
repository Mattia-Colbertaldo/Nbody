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

#include "Find_Arg.cuh"
#include "common.cuh"
#include "PhysicalForce.cuh"

#include <time.h>


#include <thrust/zip_function.h>
#include <thrust/iterator/zip_iterator.h>

constexpr float MS_PER_SEC = 1000.0f;
// constexpr int grid_size = num_parts / block_size + 1;

constexpr unsigned int BLOCK_DIM = 32;
constexpr unsigned int TILE_WIDTH = BLOCK_DIM;

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

__global__ void ResetAcc(double* ax, double* ay, double* az, const int num_parts){
  unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=num_parts) return;
  ax[i] = 0.0;
  ay[i] = 0.0;
  az[i] = 0.0;
}

class Particles{

  public:
  // Option 1
  thrust::device_vector<double3> pos;
  thrust::device_vector<double3> vel;
  thrust::device_vector<double3> acc;
  double3* dpos;
  double3* dvel;
  double3* dacc;
  // Option 2
  thrust::host_vector<double> x_h;
  thrust::host_vector<double> y_h;
  thrust::host_vector<double> z_h;
  thrust::device_vector<double> x;
  thrust::device_vector<double> y;
  thrust::device_vector<double> z;
  thrust::device_vector<double> vx;
  thrust::device_vector<double> vy;
  thrust::device_vector<double> vz;
  thrust::device_vector<double> ax;
  thrust::device_vector<double> ay;
  thrust::device_vector<double> az;
  thrust::device_vector<double> masses;
  thrust::device_vector<double> charges;
  double* dx;
  double* dy;
  double* dz;
  double* hx;
  double* hy;
  double* hz;
  double* dvx;
  double* dvy;
  double* dvz;
  double* dax;
  double* day;
  double* daz;
  double* dmasses;
  double* dcharges;
  const int num_parts;
  double size;
  AbstractForce* force;

  int th_per_block = 10;
  dim3 block_sizes;
  // block_sizes.x;
  // block_sizes.y;
  // block_sizes.z;
  // dim3 grid_sizes;
  dim3 grid_sizes;
  // grid_sizes.x;
  // grid_sizes.y;
  // grid_sizes.z;
  
  Particles(const int num_parts, AbstractForce* force) : num_parts(num_parts), force(force){
    this->size = std::sqrt(density * num_parts);
    
    thrust::default_random_engine rng;

    //Option 1
    pos.reserve(num_parts);
    vel.reserve(num_parts);
    acc.reserve(num_parts);
    // initialize_0_size(pos, rng, size);
    // initialize_11(vel, rng);
    thrust::fill(acc.begin(), acc.end(), make_double3(0.0, 0.0, 0.0));
    dpos = thrust::raw_pointer_cast(pos.data());
    dvel = thrust::raw_pointer_cast(vel.data());
    dacc = thrust::raw_pointer_cast(acc.data());
    // Option 2
    x_h.reserve(num_parts);
    y_h.reserve(num_parts);
    z_h.reserve(num_parts);
    x.reserve(num_parts);
    y.reserve(num_parts);
    z.reserve(num_parts);
    vx.reserve(num_parts);
    vy.reserve(num_parts);
    vz.reserve(num_parts);
    ax.reserve(num_parts);
    ay.reserve(num_parts);
    az.reserve(num_parts);
    masses.reserve(num_parts);
    charges.reserve(num_parts);

    // Copying to host vectors (Optional)
    
    // Block and Grid dimensions
    th_per_block = 32;
    block_sizes.x = 32;
    block_sizes.y = 32;
    // block_sizes.x = 1024;
    // block_sizes.y = 1024;
    // block_sizes.z = 0;
    // grid_sizes(ceil(num_parts/th_per_block), ceil(num_parts/th_per_block));
    grid_sizes.x = 1;
    grid_sizes.y = 1;
    // grid_sizes.x = 2147483647/1024;
    // grid_sizes.y = 65535/1024;
    // grid_sizes.z = 0;
  };

  void init(){
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
    dx = thrust::raw_pointer_cast(x.data());
    dy = thrust::raw_pointer_cast(y.data());
    dz = thrust::raw_pointer_cast(z.data());
    hx = thrust::raw_pointer_cast(x.data());
    hy = thrust::raw_pointer_cast(y.data());
    hz = thrust::raw_pointer_cast(z.data());
    dvx = thrust::raw_pointer_cast(vx.data());
    dvy = thrust::raw_pointer_cast(vy.data());
    dvz = thrust::raw_pointer_cast(vz.data());
    dax = thrust::raw_pointer_cast(ax.data());
    day = thrust::raw_pointer_cast(ay.data());
    daz = thrust::raw_pointer_cast(az.data());
    dmasses = thrust::raw_pointer_cast(masses.data());
    dcharges = thrust::raw_pointer_cast(charges.data());

    cudaDeviceSynchronize();
    ResetAcc<<<ceil(num_parts/th_per_block), th_per_block>>>(dax, day, daz, num_parts);
    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
      return;
    }
  }

  
  void save_output(std::ofstream& fsave, int step);
  void save(std::ofstream& fsave);
  void simulate_one_step();
};

void Particles::save(std::ofstream& fsave){
    
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
        fsave <<  x[i] << " " << y[i] << " " << z[i] << std::endl;
  }
}

void Particles::save_output(std::ofstream& fsave, int step){
    // TODO FIX
    // thrust::copy(x.begin(), x.end(), x_h.begin());
    // thrust::copy(y.begin(), y.end(), y_h.begin());
    // thrust::copy(z.begin(), z.end(), z_h.begin());
    save(fsave);
    if(step > 0){
        if (step%10 == 0){
        fflush(stdout);
        printf("[ %d% ]\r", (int)(step*100/nsteps));
        }
    }
}



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

  

void Particles::simulate_one_step(){
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
  move_kernel<<<ceil((double)(num_parts)/(double)1024), 1024>>>(dx, dy, dz, dvx, dvy, dvz, dax, day, daz, size, num_parts);
  cudaDeviceSynchronize();
  error = cudaGetLastError();
  if (error != cudaSuccess) {
    fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    return;
  }
  
  
  if(first) std::cout << "Moving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  if(first) first = 0;
}




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
  std::ofstream fsave(savename);
  if (savename != "") std::cout << "File created." << std::endl;

  //Find force
  std::string forcename = finder.find_string_arg("-f", "repulsive");
  if (forcename != "") std::cout << "Choosing " <<  forcename << " force..." << std::endl;
  else{
      std::string def="default";
      forcename= &def[0];
      std::cout << "Choosing default force..." << std::endl;;
  }

  AbstractForce* force= finder.find_force(forcename);

  const int num_parts = finder.find_int_arg("-n", 1000);
  std::cout << "Starting simulation with " << num_parts << " particles." <<std::endl;
  const int part_seed = finder.find_int_arg("-s", 0);
  
  // const double size = std::sqrt(density * num_parts);
  // SIZE = size;
  // std::cout << num_parts << " " << size << " " << nsteps << std::endl;
  const int num_th = finder.find_int_arg("-t", 8);

  long t = clock();

  Particles p = Particles(num_parts, force);
  p.init();
  
  cudaDeviceSynchronize();

  std::cout << "Initialization: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  // std::cout << "Moving data from GPU to CPU: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  t = clock();
  // thrust::copy(x.begin(), x.end(), x_h.begin());
  // thrust::copy(y.begin(), y.end(), y_h.begin());
  // thrust::copy(z.begin(), z.end(), z_h.begin());
  p.save(fsave);
  std::cout << "Saving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  std::cout << "Now entering the for loop." << std::endl;
  t = clock();
  long t1;
  for(int step=0; step<nsteps; step++){
    if(step == 0) t1 = clock();
    p.simulate_one_step();
    if(step == 0) std::cout << "Simulating one step: " << ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
    cudaDeviceSynchronize();
    if(step == 0) t1 = clock();
    p.save_output(fsave, step);
    if(step == 0) std::cout << "Saving: " << ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
    if(step == 0) std::cout << "One loop iteration: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  }
  std::cout << "Simulating all the steps: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  


  return 0;
}