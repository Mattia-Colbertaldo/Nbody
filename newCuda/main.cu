#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/host_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/random.h>
#include <fstream>

#include <thrust/transform.h>
#include <thrust/execution_policy.h>

#include "Find_Arg.cuh"
#include "common.cuh"
#include "PhysicalForce.cuh"



#include <cassert>
#include <iostream>

double SIZE;


void initialize(thrust::universal_vector<double>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(0., SIZE);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng);

  keys = h_keys;
}

void initialize_1(thrust::universal_vector<double>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng);

  keys = h_keys;
}

void initialize_2(thrust::universal_vector<double>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng)*1e-19;

  keys = h_keys;
}

void save(std::ofstream& fsave, thrust::universal_vector<double>& x, thrust::universal_vector<double>& y, thrust::universal_vector<double>& z, int num_parts){
    
  static bool first = true;

  if (first) {
      fsave << num_parts << " " << SIZE << " " << nsteps << "\n";
      first = false;
  }

  for(size_t i = 0; i < num_parts; i++){
        fsave <<  x[i] << " " << y[i] << " " << z[i] << std::endl;
  }
}

void save_output(std::ofstream& fsave, thrust::universal_vector<double>& x, thrust::universal_vector<double>& y, thrust::universal_vector<double>& z,
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
                __host__ __device__
                    double operator()(const double& v, const double& x) const {
                        return dt * v + x;
                    }
            };


            void saxpy_fast(thrust::universal_vector<double>& vx, thrust::universal_vector<double>& x)
            {
              // x <- dt * vx + x
              thrust::transform(thrust::device, vx.begin(), vx.end(), x.begin(), x.begin(), saxpy_functor());
            }

            void move(thrust::universal_vector<double>& x, thrust::universal_vector<double>& y, thrust::universal_vector<double>& z,
                      thrust::universal_vector<double>& vx, thrust::universal_vector<double>& vy, thrust::universal_vector<double>& vz,
                      thrust::universal_vector<double>& ax, thrust::universal_vector<double>& ay, thrust::universal_vector<double>& az){
                    saxpy_fast(ax, vx);
                    saxpy_fast(ay, vy);
                    saxpy_fast(az, vz);
                    saxpy_fast(vx, x);
                    saxpy_fast(vy, y);
                    saxpy_fast(vz, z);
                  }

              // **************************** APPLYING FORCE **************************** //

            void apply_force(thrust::universal_vector<double>& x, thrust::universal_vector<double>& y, thrust::universal_vector<double>& z,
                             thrust::universal_vector<double>& vx, thrust::universal_vector<double>& vy, thrust::universal_vector<double>& vz,
                             thrust::universal_vector<double>& ax, thrust::universal_vector<double>& ay, thrust::universal_vector<double>& az){

                    //TODO: MODIFY, JUST FOR TESTING
                    initialize_1(ax);
                    initialize_1(ay);
                    initialize_1(az);
            }




void simulate_one_step(thrust::universal_vector<double>& x, thrust::universal_vector<double>& y, thrust::universal_vector<double>& z,
                  thrust::universal_vector<double>& vx, thrust::universal_vector<double>& vy, thrust::universal_vector<double>& vz,
                  thrust::universal_vector<double>& ax, thrust::universal_vector<double>& ay, thrust::universal_vector<double>& az,
                  AbstractForce& force){
                  thrust::fill(ax.begin(), ax.end(), 0.);
                  thrust::fill(ay.begin(), ay.end(), 0.);
                  thrust::fill(az.begin(), az.end(), 0.);
                  cudaDeviceSynchronize();
                  apply_force(x, y, z, vx, vy, vz, ax, ay, az);
                  cudaDeviceSynchronize();
                  move(x, y, z, vx, vy, vz, ax, ay, az);
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

  thrust::universal_vector<double> x(num_parts);
  thrust::universal_vector<double> y(num_parts);
  thrust::universal_vector<double> z(num_parts);
  thrust::universal_vector<double> vx(num_parts);
  thrust::universal_vector<double> vy(num_parts);
  thrust::universal_vector<double> vz(num_parts);
  thrust::universal_vector<double> ax(num_parts);
  thrust::universal_vector<double> ay(num_parts);
  thrust::universal_vector<double> az(num_parts);
  thrust::universal_vector<double> masses(num_parts);
  thrust::universal_vector<double> charges(num_parts);
  initialize(x);
  initialize(y);
  initialize(z);
  initialize_2(vx);
  initialize_2(vy);
  initialize_2(vz);
  thrust::fill(ax.begin(), ax.end(), 0.);
  thrust::fill(ay.begin(), ay.end(), 0.);
  thrust::fill(az.begin(), az.end(), 0.);
  initialize_1(masses);
  initialize_2(charges);

  double* px = thrust::raw_pointer_cast(x.data());
  double* py = thrust::raw_pointer_cast(y.data());
  double* pz = thrust::raw_pointer_cast(z.data());
  double* pvx = thrust::raw_pointer_cast(vx.data());
  double* pvy = thrust::raw_pointer_cast(vy.data());
  double* pvz = thrust::raw_pointer_cast(vz.data());
  double* pax = thrust::raw_pointer_cast(ax.data());
  double* pay = thrust::raw_pointer_cast(ay.data());
  double* paz = thrust::raw_pointer_cast(az.data());
  double* pmasses = thrust::raw_pointer_cast(masses.data());
  double* pcharges = thrust::raw_pointer_cast(charges.data());

  save(fsave, x, y, z, num_parts);
  for(int step=0; step<nsteps; step++){
    simulate_one_step(x, y, z, vx, vy, vz, ax, ay, az, *force);
    cudaDeviceSynchronize();
    
    save_output(fsave, x, y, z, step, num_parts);
  }
  


  return 0;
}
