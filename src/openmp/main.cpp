#include "common.h"
#include "Find_Arg.hpp"
#include "PhysicalForce.hpp"
#include "Output.hpp"
#include "Particle.hpp"
#include "Simulation.hpp"

#include<cmath>
#include <chrono>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <thread>
#include <omp.h>
#include <stdexcept>


using namespace common_h;


int main(int argc, char** argv) {
    // Parse Args
    Find_Arg finder= Find_Arg(argc, argv);
    finder.find_int_arg("-h", 0);

    // Open Output File
    std::string savename = finder.find_string_option("-o", "out.txt");
    if (savename != "") std::cout << "Output: " << savename << "..." << std::endl;
    std::ofstream fsave(savename);

    //Find force
    std::string forcename = finder.find_string_option("-f", "repulsive");

    std::unique_ptr<AbstractForce> force= finder.find_force(forcename);

    //find collision type
    int collision = finder.find_int_arg("-c", 0);
    std::cout << "Collision: " <<  collision << std::endl;



    // Initialize Particles
    const int num_parts = finder.find_int_arg("-n", 1000);
    const int part_seed = finder.find_int_arg("-s", 0);
   
    const double size = std::sqrt(density * num_parts);
    std::cout << num_parts << " " << size << " " << nsteps << std::endl;
    const int num_th = finder.find_int_arg("-t", 8);

    
    Simulation simulation = Simulation(num_parts, collision);
    std::cout << "Initialization: ";
    simulation.init_particles(size, part_seed);
    
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    auto init_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff_1 = init_time - start_time;
    double seconds_1 = diff_1.count();
    std::cout << seconds_1 << " seconds" << std::endl;

Output output = Output(savename);
output.save(simulation.parts , size);
#ifdef _OPENMP
std::cout << "Available threads: " << std::thread::hardware_concurrency() << "\nRunning "
          << num_th << " thread" << (num_th>1? "s." : ".") <<std::endl;
#pragma omp parallel default(shared) num_threads(num_th)
#endif
    {
        for (int step = 0; step < nsteps; ++step) {
            
            simulation.simulate_one_step(force, num_parts,size);

            // Save state if necessary
            #ifdef _OPENMP
            #pragma omp master
            #endif
            {
                output.save_output(simulation.parts , step, size);
            }
            
        }
        
        
    }

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    // Finalize
    std::cout << "Simulation Time = " << seconds << " seconds for " << num_parts <<
     " particles and " << nsteps << " steps.\n";
    fsave.close();

}
