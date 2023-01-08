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

constexpr unsigned int TILE_WIDTH = 10;
constexpr unsigned int BLOCK_DIM = 10;


#include <cassert>
#include <iostream>

#define BLOCK_SIZE = 32;




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

  void apply_force();
  
  void save_output(std::ofstream& fsave, int step);
  void save(std::ofstream& fsave);
  void simulate_one_step();
  void move();

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

            // **************************** MOVING **************************** //

struct saxpy_functor
{
saxpy_functor(const double size) : size(size){};
double size;
__device__
    double operator()(double& v, double& x) const {
        double result = dt * v + x;
        // Bounce from walls
        if(result < 0 || result > size){
          result = result < 0 ? -result : 2*size - result;
          v = -v;
        }
        return result;
    }
};

struct saxpy_functor2
{
saxpy_functor2(const double size) : size(size){};
double size;
__device__
    double operator()(double& a, double& v) const {
        return dt * a + v;
    }
};


void saxpy_fast(thrust::device_vector<double>& vx, thrust::device_vector<double>& x, const double size)
{
// x <- dt * vx + x
thrust::transform(thrust::device, vx.begin(), vx.end(), x.begin(), x.begin(), saxpy_functor(size));
}

void saxpy_fast2(thrust::device_vector<double>& ax, thrust::device_vector<double>& vx, const double size)
{
// vx <- dt * ax + vx
thrust::transform(thrust::device, ax.begin(), ax.end(), vx.begin(), vx.begin(), saxpy_functor2(size));
}

void Particles::move(){
    saxpy_fast2(ax, vx, size);
    saxpy_fast2(ay, vy, size);
    saxpy_fast2(az, vz, size);
    cudaDeviceSynchronize();
    saxpy_fast(vx, x, size);
    saxpy_fast(vy, y, size);
    saxpy_fast(vz, z, size);
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
    if (dx[i] < 0 || dx[i] > size) {
        dx[i] = (dx[i] < 0 ? -dx[i] : 2 * size - dx[i]);
        dvx[i] = -dvx[i];
    }

    if (dy[i] < 0 || dy[i] > size) {
        dy[i] = (dy[i] < 0 ? -dy[i] : 2 * size - dy[i]);
        dvy[i] = -dvy[i];
    }

    if (dz[i] < 0 || dz[i] > size) {
        dz[i] = (dz[i] < 0 ? -dz[i] : 2 * size - dz[i]);
        dvz[i] = -dvz[i];
    }
}

// **************************** APPLYING FORCE **************************** //

struct arbitrary_functor
{
arbitrary_functor(double& x_i, double& y_i, double& z_i, double& ax_i, double& ay_i, double& az_i) :
  x_i(x_i), y_i(y_i), z_i(z_i), ax_i(ax_i), ay_i(ay_i), az_i(az_i){};
double x_i, y_i, z_i, ax_i, ay_i, az_i;
__device__
void operator()(const double& x_j, const double& y_j, const double& z_j, const double& m_j)
{   double dx = x_j - x_i;
    double dy = y_j - y_i;
    double dz = z_j - z_i;
    double r2_j = dx * dx + dy * dy + dz * dz;
    if (r2_j > cutoff * cutoff){
      ax_i = 0.0;
      ay_i = 0.0;
      az_i = 0.0;
    }
    else{
      r2_j = fmax(r2_j, min_r * min_r);
      double coef = G * ( m_j / r2_j );
      ax_i = coef*dx;
      ay_i = coef*dx;
      az_i = coef*dx;
      // d += a + b * c;
    }
    
}
};

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

__global__ void force_kernel(double* x, double* y, double* z,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts){
    // __shared__ double tile[250];
    int thx = threadIdx.x;
    int bx = blockIdx.x;
    int index = bx * blockDim.x + thx;
    __shared__ double tilex[TILE_WIDTH];
    __shared__ double tiley[TILE_WIDTH];
    __shared__ double tilez[TILE_WIDTH];
    __shared__ double tile_masses[TILE_WIDTH];
    __shared__ double tile_charges[TILE_WIDTH];
    
    // int index = blockIdx.x * 250 + thx;

    // unsigned int i = blockDim.y * blockIdx.y + threadIdx.y;
    
    // I do calculations only if x index is a successor of y index

    // Loop over all particles with tiles
    for(unsigned int i=0; i< (num_parts-1)/(TILE_WIDTH)+1 ; i++){ 
      // Collaborative loading
      if(thx < 0 || thx >= TILE_WIDTH) return;
      if(index < num_parts &&  (i* TILE_WIDTH+thx) < num_parts) {
        // tilex[thx] = x[ i*TILE_WIDTH+thx ];
        tilex[thx] = x[ index ];
        tiley[thx] = y[ index ];
        tilez[thx] = z[ index ];
        tile_masses[thx] = masses[ index ];
        tile_charges[thx] = charges[ index ];
      }
      else
      {
        tilex[thx] = 0.0;
        tiley[thx] = 0.0;
        tilez[thx] = 0.0;
        tile_masses[thx] = 0.0;
        tile_charges[thx] = 0.0;
      }

      __syncthreads();
      // here every thread calculates the force between his particle and the particles in the tile
      if(index < num_parts){
        // loop inside the tile
        for(unsigned int j=0; j <TILE_WIDTH; j++){
          double dx = x[index] - tilex[j];
          double dy = y[index] - tiley[j];
          double dz = z[index] - tilez[j];
          double r2 = dx * dx + dy * dy + dz * dz;
          if (r2 > cutoff * cutoff) {return;}
          else{
            r2 = fmax(r2, min_r * min_r);
            double coef = G * tile_masses[j] / r2;

            // 1
            atomicAdd((double*)(ax + index), (double)coef*dx);
            // ax[index] += coef*dx;
            atomicAdd((double*)(ay + index), (double)coef*dy);
            // ay[index] += coef*dy;
            ay[index]= ay[index] + coef*dy;
            atomicAdd((double*)(az + index), (double)coef*dz);
            // az[index] += coef*dz;
            az[index]= az[index] + coef*dz;
          
          }
        }
      }
      __syncthreads();
  }

    
}

__global__ void kernel_test_force(double* x, double* y, double* z, double* vx, double* vy, double* vz,
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

void Particles::simulate_one_step(){
  static bool first = 1;
  long t;
  if(first) t = clock();
  
  // long t = clock();
  // apply_force(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts, sum_ax, sum_ay, sum_az);
  long t1;
  if(first) t1 = clock();
  th_per_block = fmin(32, num_parts);
  block_sizes.x = block_sizes.y = 32;
  block_sizes.z = 1;
  grid_sizes.x = ceil(((double)num_parts)/((double)(block_sizes.x)));
  grid_sizes.y = ceil(((double)num_parts)/((double)(block_sizes.y)));
  grid_sizes.z = 1;
  if(first) std::cout << "GRID SIZE: " << grid_sizes.x << std::endl;
  ResetAcc<<<ceil((double)(num_parts)/(double)1024), 1024>>>(dax, day, daz, num_parts);
  kernel_test_force<<<grid_sizes, block_sizes>>>(dx, dy, dz, dvx, dvy, dvz, dax, day, daz, dmasses, dcharges, num_parts);
  // force_kernel<<<ceil((double)(num_parts)/(double)1024), 1024>>>(dx, dy, dz, dax, day, daz, dmasses, dcharges, num_parts);
  if(first) std::cout << "Applying force: Kernel Loop: " << ((clock() - t1)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  
  cudaDeviceSynchronize();
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    return;
  }
  
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
  // move();
  // cudaDeviceSynchronize();
  // error = cudaGetLastError();
  // if (error != cudaSuccess) {
  //   fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
  //   return;
  // }
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
