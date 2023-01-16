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

// ==============
// Main Function
// ==============



int main(int argc, char** argv) {
    // Parse Args
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

    std::unique_ptr<AbstractForce> force= finder.find_force(forcename);

    //find collision type
    int collision = finder.find_int_arg("-c", 0);
    std::cout << "Choosing " <<  collision << " collision type..." << std::endl;


    

    
    /*
    const std::unordered_map<std::string, std::shared_ptr<AbstractForce>> fmap =
    {
        {"default", std::make_shared<RepulsiveForce>() },
        {"repulsive", std::make_shared<RepulsiveForce>() },
        {"gravitational", std::make_shared<GravitationalForce>() },
        {"assist", std::make_shared<GravitationalAssistForce>() },
        {"proton", std::make_shared<ProtonForce>() },
        {"coulomb", std::make_shared<CoulombForce>() },
    };
    std::shared_ptr<AbstractForce> force;

    try {
        force = fmap.at(forcename);      // vector::at throws an out-of-range
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "Wrong Force chosen "<< '\n';
        force =  fmap.at("repulsive");
    }

    */



    
    

    // Initialize Particles
    const int num_parts = finder.find_int_arg("-n", 1000);
    const int part_seed = finder.find_int_arg("-s", 0);
   
    const double size = std::sqrt(density * num_parts);
    std::cout << num_parts << " " << size << " " << nsteps << std::endl;
    const int num_th = finder.find_int_arg("-t", 8);

    
    Simulation simulation = Simulation(num_parts, collision);
    std::cout << "Trying to init particles..." << std::endl;
    simulation.init_particles(size, part_seed);
    
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    std::cout << "Trying to init simulation..." << std::endl;
    std::cout << "Init simulation ended." << std::endl;

    auto init_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff_1 = init_time - start_time;
    double seconds_1 = diff_1.count();
    std::cout << "initialization Time = " << seconds_1 << " seconds\n";

Output output = Output();
output.save(fsave, simulation.parts , size, nsteps);
#ifdef _OPENMP
std::cout << "Available threads: " << std::thread::hardware_concurrency() << "\nRunning "
          << num_th << " thread(s)." <<std::endl;
#pragma omp parallel default(shared) num_threads(num_th)
#endif
    {
        //for nel tempo: non parallelizzare
        for (int step = 0; step < nsteps; ++step) {
            
            simulation.simulate_one_step(force, num_parts,size);

            // Save state if necessary
            #ifdef _OPENMP
            #pragma omp master
            #endif
            {
                output.save_output(fsave, savefreq, simulation.parts , step, nsteps, size);
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
