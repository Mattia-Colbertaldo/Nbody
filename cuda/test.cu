#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

#include "common.cuh"
#include "Find_Arg.cuh"
#include "PhysicalForce.cuh"
#include "Output.cuh"
#include "Particles.cuh"

#include<cmath>
#include <chrono>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <thread>
#include <omp.h>
#include <stdexcept>


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

#define MATRIX_SIZE 1000*1000
#define CPU_VECTOR_SIZE 1024
#define MASK_WIDTH 221
#define TILE_WIDTH 221

double SIZE;


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

__global__ void simulate_one_step(Particles particles){
    particles.simulate_one_step();
}

__global__ void save(Output output, const char* charname, const int num_parts, const Particles& particles, const double size, const int& nsteps){
    output.save(charname, num_parts, particles , size, nsteps);
}

__global__ void save_output(Output output, const char* charname, const int savefreq, const int num_parts, const Particles& particles , const int& step,  const int& nsteps, const double & size){
    output.save_output(charname, savefreq, num_parts, particles, step, nsteps, size);
}

__global__ void gpuPrint(int* i){
    printf("GPU:  %d\n", *i);
}

__global__ void gpuPrint(double* i){
    printf("GPU:  %f\n", *i);
}

void initialize(thrust::device_vector<double3>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(0., SIZE);

  thrust::host_vector<double3> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = make_double3(dist(rng), dist(rng), dist(rng));

  keys = h_keys;
}

void initialize_1(thrust::device_vector<double3>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double3> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = make_double3(dist(rng), dist(rng), dist(rng));

  keys = h_keys;
}

void initialize_2(thrust::device_vector<double>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng);

  keys = h_keys;
}

void initialize_3(thrust::device_vector<double>& keys)
{
  thrust::default_random_engine rng;
  thrust::uniform_real_distribution<double> dist(-1.0, 1.0);

  thrust::host_vector<double> h_keys(keys.size());

  for(size_t i = 0; i < h_keys.size(); i++)
    h_keys[i] = dist(rng)*1e-19;

  keys = h_keys;
}




// ==============
// Main Function
// ==============



int main(int argc, char** argv) {

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
    const char * charname = savename.c_str();
    if (savename != "") std::cout << "Creating file " << savename << "..." << std::endl;
    // std::ofstream fsave(savename);
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

    int num_parts = finder.find_int_arg("-n", 1000);
    const int part_seed = finder.find_int_arg("-s", 0);
   
    double size = std::sqrt(density * num_parts);
    SIZE = size;
    std::cout << num_parts << " " << size << " " << nsteps << std::endl;
    const int num_th = finder.find_int_arg("-t", 8);






    // **************************INITIALIZATION ENDED************************** //
    
    thrust::device_vector<double3> positions(num_parts);
    thrust::device_vector<double3> velocities(num_parts);
    thrust::device_vector<double3> accelerations(num_parts);
    thrust::device_vector<double> masses(num_parts);
    thrust::device_vector<double> charges(num_parts);
    initialize(positions);
    initialize_1(velocities);
    initialize_2(masses);
    initialize_3(charges);
    thrust::fill(accelerations.begin(), accelerations.end(), make_double3(0., 0., 0.));

    //std::cout << positions[0] << std::endl;

    return;

    // thrust::device_ptr<Particles> p = thrust::device_malloc<Particles>(1);

    
    
    
    // // device_ptr can be converted to a "raw" pointer for use in other APIs and kernels, etc.
    // Particles * raw_ptr = thrust::raw_pointer_cast(p);

    // // note: raw_ptr cannot necessarily be accessed by the host!

    // // conversely, raw pointers can be wrapped
    // thrust::device_ptr<Particles> wrapped_ptr = thrust::device_pointer_cast(raw_ptr);

    // raw_ptr = new Particles(num_parts, size);

    // std::cout << (*raw_ptr).x[0] << std::endl;

    // std::cout << (*raw_ptr).x[0] << std::endl;

    


    // std::cout << "Initializing..." << std::endl;
    // particles.Initialize();
    // cudaDeviceSynchronize();
    // std::cout << "Setting the force..." << std::endl;
    // particles.setForce(force);
    // cudaDeviceSynchronize();
    // std::cout << "Force set." << std::endl;

    // cudaEvent_t start, stop;
    // float elapsedTime;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // cudaEventRecord(start, 0);
    
    // std::cout << "Redirecting stdout to " << charname << "..." <<std::endl;
    // // fclose (stdout);
    // // freopen (charname, "w", stdout);
    // Output output = Output();
    // std::cout << "Saving first output..." << std::endl;
    // save<<<1,1>>>(output, charname, num_parts, particles, size, nsteps);
    // cudaDeviceSynchronize();
    // std::cout << "Saved." << std::endl;

    // //for nel tempo: non parallelizzare
    // for (int step = 0; step < nsteps; ++step) {
    //     std::cout << "Simulating step " << step << std::endl;
    //     // TODO: CHANGE INTO KERNEL
    //     simulate_one_step<<<1, 1>>>(particles);
    //     std::cout << "Step " << step << " simulation started, waiting for gpu to finish..." << std::endl;
    //     cudaDeviceSynchronize();
    //     std::cout << "GPU synchronized." << std::endl;

    //     // Save state if necessary
    //     {
    //         std::cout << "Saving output..." << std::endl;
    //         save_output<<<1,1>>>(output, charname, num_parts, savefreq, particles, step, nsteps, size);
    //         cudaDeviceSynchronize();
    //         std::cout << "Saved." << std::endl;
    //     }
        
    // }
        
    
    // cudaThreadSynchronize();

    // // time counting terminate
    
    // cudaEventRecord(stop, 0);
    // cudaEventSynchronize(stop);

    // // compute time elapsed on GPU computing
    // cudaEventRecord(stop,0);
    // cudaEventSynchronize(stop);
    // cudaEventElapsedTime(&elapsedTime, start, stop);

    // // Finalize
    // std::cout << "Simulation Time = " << elapsedTime/1000 << " seconds for " << num_parts <<
    //  " particles and " << nsteps << " steps.\n";
    // fflush (stdout);
    // fclose (stdout);

    // free memory
    // cudaFree(&simulation);
    // cudaFree(&force);
    // cudaFree(&num_parts);
    // cudaFree(&size);

}
