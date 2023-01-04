#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

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



#include <cassert>
#include <iostream>

double SIZE;


void initialize(thrust::device_vector<double>& keys, thrust::default_random_engine& rng, const double& size)
{
  thrust::uniform_real_distribution<double> dist(0., SIZE);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng);

  keys = h_keys;
}

void initialize_1(thrust::device_vector<double>& keys, thrust::default_random_engine& rng)
{
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng);

  keys = h_keys;
}

void initialize_2(thrust::device_vector<double>& keys, thrust::default_random_engine& rng)
{
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng)*1e-19;

  keys = h_keys;
}

void save(std::ofstream& fsave, thrust::device_vector<double>& x, thrust::device_vector<double>& y, thrust::device_vector<double>& z, int num_parts){
    
  static bool first = true;

  if (first) {
      fsave << num_parts << " " << SIZE << " " << nsteps << "\n";
      first = false;
  }

  for(size_t i = 0; i < num_parts; i++){
        fsave <<  x[i] << " " << y[i] << " " << z[i] << std::endl;
  }
}

void save_output(std::ofstream& fsave, thrust::device_vector<double>& x, thrust::device_vector<double>& y, thrust::device_vector<double>& z,
                 int step, int num_parts){
    save(fsave, x, y, z, num_parts);
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
                saxpy_functor(const double& size) : size(size){};
                double size;
                __host__ __device__
                    double operator()(double& v, const double& x) const {
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
                saxpy_functor2(const double& size) : size(size){};
                double size;
                __host__ __device__
                    double operator()(double& a, const double& v) const {
                        return dt * a + v;
                    }
            };


            void saxpy_fast(thrust::device_vector<double>& vx, thrust::device_vector<double>& x, const double& size)
            {
              // x <- dt * vx + x
              thrust::transform(thrust::device, vx.begin(), vx.end(), x.begin(), x.begin(), saxpy_functor(size));
            }

            void saxpy_fast2(thrust::device_vector<double>& ax, thrust::device_vector<double>& vx, const double& size)
            {
              // x <- dt * vx + x
              thrust::transform(thrust::device, ax.begin(), ax.end(), vx.begin(), vx.begin(), saxpy_functor2(size));
            }

            void move(thrust::device_vector<double>& x, thrust::device_vector<double>& y, thrust::device_vector<double>& z,
                      thrust::device_vector<double>& vx, thrust::device_vector<double>& vy, thrust::device_vector<double>& vz,
                      thrust::device_vector<double>& ax, thrust::device_vector<double>& ay, thrust::device_vector<double>& az, const double& size){
                    saxpy_fast2(ax, vx, size);
                    saxpy_fast2(ay, vy, size);
                    saxpy_fast2(az, vz, size);
                    cudaDeviceSynchronize();
                    saxpy_fast(vx, x, size);
                    saxpy_fast(vy, y, size);
                    saxpy_fast(vz, z, size);
                  }

              // **************************** APPLYING FORCE **************************** //
            
            struct arbitrary_functor
            {
                arbitrary_functor(double& x_i, double& y_i, double& z_i) : x_i(x_i), y_i(y_i), z_i(z_i){};
                double x_i, y_i, z_i;
                __host__ __device__
                void operator()(const double& x_j, const double& y_j, const double& z_j, const double& m_j, double& ax_j, double& ay_j, double& az_j)
                {   double dx = x_j - x_i;
                    double dy = y_j - y_i;
                    double dz = z_j - z_i;
                    double r2_j = dx * dx + dy * dy + dz * dz;
                    if (r2_j > cutoff * cutoff){
                      ax_j = 0.;
                      ay_j = 0.;
                      az_j = 0.;
                    }
                    else{
                      r2_j = fmax(r2_j, min_r * min_r);
                      double coef = G * ( m_j / r2_j );
                      ax_j = coef*dx;
                      ay_j = coef*dx;
                      az_j = coef*dx;
                      // d += a + b * c;
                    }
                    
                }
            };


            void apply_force(thrust::device_vector<double>& x, thrust::device_vector<double>& y, thrust::device_vector<double>& z,
                             thrust::device_vector<double>& vx, thrust::device_vector<double>& vy, thrust::device_vector<double>& vz,
                             thrust::device_vector<double>& ax, thrust::device_vector<double>& ay, thrust::device_vector<double>& az,
                             thrust::device_vector<double>& masses, thrust::device_vector<double>& charges, const int& num_parts){
                    
                    thrust::device_vector<double> sum_ax(num_parts);
                    thrust::device_vector<double> sum_ay(num_parts);
                    thrust::device_vector<double> sum_az(num_parts);
                    // std::cout << "   Apply_Force --> Starting for loop" << std::endl;
                    for(int i=0; i<num_parts; i++){
                      double x_i = x[i];
                      double y_i = y[i];
                      double z_i = z[i];
                      thrust::for_each(thrust::device, thrust::make_zip_iterator(thrust::make_tuple(x.begin(), y.begin(), z.begin(), masses.begin(), sum_ax.begin(), sum_ay.begin(), sum_az.begin())),
                                        thrust::make_zip_iterator(thrust::make_tuple(x.end(), y.end(), z.end(), masses.end(), sum_ax.end(), sum_ay.end(), sum_az.end())),
                                        thrust::make_zip_function(arbitrary_functor(x_i, y_i, z_i)));
                      cudaDeviceSynchronize();
                      // std::cout << "   Apply_Force --> Reducing" << std::endl;
                      ax[i] = thrust::reduce(thrust::device, sum_ax.begin(), sum_ax.end(), 0.);
                      ay[i] = thrust::reduce(thrust::device, sum_ay.begin(), sum_ay.end(), 0.);
                      az[i] = thrust::reduce(thrust::device, sum_az.begin(), sum_az.end(), 0.);
                    }

                    
            }




void simulate_one_step(thrust::device_vector<double>& x, thrust::device_vector<double>& y, thrust::device_vector<double>& z,
                  thrust::device_vector<double>& vx, thrust::device_vector<double>& vy, thrust::device_vector<double>& vz,
                  thrust::device_vector<double>& ax, thrust::device_vector<double>& ay, thrust::device_vector<double>& az,
                  thrust::device_vector<double>& masses, thrust::device_vector<double>& charges, AbstractForce& force, const int& num_parts, const double& size){
                  thrust::fill(ax.begin(), ax.end(), 0.);
                  thrust::fill(ay.begin(), ay.end(), 0.);
                  thrust::fill(az.begin(), az.end(), 0.);
                  // std::cout << "Accelerations filled with 0s" << std::endl;
                  cudaDeviceSynchronize();
                  apply_force(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, num_parts);
                  // std::cout << "Force applied" << std::endl;
                  cudaDeviceSynchronize();
                  move(x, y, z, vx, vy, vz, ax, ay, az, size);
                  // std::cout << "Particles moved" << std::endl;
}


int main(int argc, char** argv)
{

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
  if (forcename != "") std::cout << "Choosing non default force " <<  forcename << "..." << std::endl;
  else{
      std::string def="default";
      forcename= &def[0];
      std::cout << "Choosing default force..." << std::endl;;
  }

  AbstractForce* force= finder.find_force(forcename);

  const int num_parts = finder.find_int_arg("-n", 1000);
  const int part_seed = finder.find_int_arg("-s", 0);
  
  const double size = std::sqrt(density * num_parts);
  SIZE = size;
  std::cout << num_parts << " " << size << " " << nsteps << std::endl;
  const int num_th = finder.find_int_arg("-t", 8);

  long t = clock();

  thrust::device_vector<double> x_h(num_parts);
  thrust::device_vector<double> y_h(num_parts);
  thrust::device_vector<double> z_h(num_parts);

  thrust::device_vector<double> x(num_parts);
  thrust::device_vector<double> y(num_parts);
  thrust::device_vector<double> z(num_parts);
  thrust::device_vector<double> vx(num_parts);
  thrust::device_vector<double> vy(num_parts);
  thrust::device_vector<double> vz(num_parts);
  thrust::device_vector<double> ax(num_parts);
  thrust::device_vector<double> ay(num_parts);
  thrust::device_vector<double> az(num_parts);
  thrust::device_vector<double> masses(num_parts);
  thrust::device_vector<double> charges(num_parts);
  thrust::default_random_engine rng;
  initialize(x, rng, size);
  initialize(y, rng, size);
  initialize(z, rng, size);
  initialize_2(vx, rng);
  initialize_2(vy, rng);
  initialize_2(vz, rng);
  thrust::fill(ax.begin(), ax.end(), 0.);
  thrust::fill(ay.begin(), ay.end(), 0.);
  thrust::fill(az.begin(), az.end(), 0.);
  initialize_1(masses, rng);
  initialize_2(charges, rng);
  cudaDeviceSynchronize();

  std::cout << "Initialization: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  
  t = clock();
  thrust::copy(x.begin(), x.end(), x_h.begin());
  thrust::copy(y.begin(), y.end(), y_h.begin());
  thrust::copy(z.begin(), z.end(), z_h.begin());
  std::cout << "Moving data from GPU to CPU: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  t = clock();
  save(fsave, x_h, y_h, z_h, num_parts);
  std::cout << "Saving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  t = clock();
  for(int step=0; step<nsteps; step++){
    simulate_one_step(x, y, z, vx, vy, vz, ax, ay, az, masses, charges, *force, num_parts, size);
    if(step == 1) std::cout << "Simulating one step: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
    cudaDeviceSynchronize();
    thrust::copy(x.begin(), x.end(), x_h.begin());
    thrust::copy(y.begin(), y.end(), y_h.begin());
    thrust::copy(z.begin(), z.end(), z_h.begin());
    save_output(fsave, x_h, y_h, z_h, step, num_parts);
  }
  std::cout << "Simulating all the steps: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
  


  return 0;
}
