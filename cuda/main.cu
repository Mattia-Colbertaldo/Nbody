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


int main(int argc, char** argv)
{

  DisplayHeader();

  Find_Arg finder= Find_Arg(argc, argv);
  finder.find_string_option("-h", 0);

  // Open Output File
  std::string savename = finder.find_string_option("-o", "out.txt");
  if (savename != "") std::cout << "Output: " << savename << std::endl;
  
  //Find force
  std::string forcename = finder.find_string_option("-f", "repulsive");

  std::shared_ptr<AbstractForce> force =finder.find_force(forcename) ;

  //find collision type
  int collision = finder.find_string_option("-c", 0);
  std::cout << "Collision: " <<  collision << std::endl;


  const int num_parts = finder.find_string_option("-n", 1000);
  std::cout << "Starting simulation with " << num_parts << " particles." <<std::endl;
  const int part_seed = finder.find_string_option("-s", 0);
  
  std::unique_ptr<AllParticles> all_particles = std::make_unique<AllParticles>(num_parts, force);
  
  const double size = std::sqrt(density * num_parts);
  const int num_th = finder.find_string_option("-t", 8);
  long t = clock();
  Simulation s = Simulation(all_particles, collision );
  std::cout << "Initialization: ";
  s.init_particles(size, part_seed);
  cudaDeviceSynchronize();
  

  std::cout << ((clock() - t)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;
  t = clock();
  Output output= Output(num_parts, savename);
  output.save( s.parts, size);
  std::cout << "Saving: " << ((clock() - t)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;
  std::cout << "Now entering the for loop." << std::endl;
  t = clock();
  long t1;
  long savetime = 0.;
  for(int step=0; step<nsteps; step++){
    if(step == nsteps/2) t1 = clock();
    s.simulate_one_step(force, num_parts, size);
    cudaDeviceSynchronize();
    if(step == nsteps/2) std::cout << "Simulating one step: " << ((clock() - t1)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;
    if(step == nsteps/2) t1 = clock();
    output.save_output( s.parts , step, size);
    if(step == nsteps/2) std::cout << "Saving to device buffer: " << ((clock() - t1)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;
    else if(step == nsteps - 1){
      savetime = ((clock() - t1)*1000)/CLOCKS_PER_SEC;
      std::cout << "Final Saving: " << savetime << std::endl;
    }
    if(step == nsteps/2) std::cout << "One loop iteration: " << savetime << " ms" << std::endl;
  }
  long all_steps = ((clock() - t)*1000)/CLOCKS_PER_SEC;
  std::cout << std::endl << "***************************************" << std::endl;
  std::cout << " + Calculations: " << all_steps-savetime << " ms [ " << ((all_steps-savetime)*100/all_steps) << " % ]" << std::endl;
  std::cout << " + Saving: " << savetime << " ms [ " << (savetime*100/all_steps) << " % ]" << std::endl;
  std::cout << " = Simulating all the steps: " << all_steps << " ms" << std::endl;
  std::cout << "***************************************" << std::endl << std::endl;;
  


  return 0;
}
