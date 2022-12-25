#include "common.h"
#include "Force.hpp"
#include "particle.hpp"
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



// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream& fsave, std::vector<particle>& parts, double size) {
    int num_parts = parts.size();
    static bool first = true;

    if (first) {
        fsave << num_parts << " " << size << " " << nsteps << "\n";
        first = false;
    }

    for (int i = 0; i < num_parts; ++i) {
        fsave << parts[i].x << " " << parts[i].y << " " << parts[i].z << "\n";
    }

    // fsave << std::endl;
}

// Particle Initialization
void init_particles(std::vector<particle>& parts, double size, int part_seed) {
    /*
        input :  1. parts     : vettore di particelle
                 2. masses    : vettore delle masse
                 3. size      : dimensione della particella
                 4. part_seed : seme randomico
    */

    int num_parts = parts.size();

    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());

    int sx = (int)ceil(sqrt((double)num_parts));
    int sy = (num_parts + sx - 1) / sx;
    int sz = sy;

    std::vector<int> shuffle(num_parts);
    for (int i = 0; i < shuffle.size(); ++i) {
        shuffle[i] = i;
    }

    for (int i = 0; i < num_parts; ++i) {
        
        // Make sure particles are not spatially sorted
        std::uniform_int_distribution<int> rand_int(0, num_parts - i - 1);
        int j = rand_int(gen);
        int k = shuffle[j];
        shuffle[j] = shuffle[num_parts - i - 1];

        //
        const double x = size * (1. + (k % sx)) / (1 + sx);
        const double y = size * (1. + (k / sx)) / (1 + sy);
        const double z = x*y;

        //
        std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
        const double vx = rand_real(gen);
        const double vy = rand_real(gen);
        const double vz = rand_real(gen);

        //
        std::uniform_real_distribution<float> rand_mass(0.001, 0.1);
        const float m = rand_mass(gen);
        parts[i]=particle(x, y, z, vx, vy, vz, m);
    }

}


void simulate_one_step(std::vector<particle>& parts, const std::shared_ptr<AbstractForce> force, int num_parts, double size) {
    // Compute Forces
    //int num_parts = parts.size();
    #ifdef _OPENMP
	#pragma omp for schedule(dynamic)
    #endif
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = parts[i].az = 0.;
        for (int j = 0; j < num_parts; ++j) {
            parts[i].apply_force(parts[j], force);
        }
    }
    #ifdef _OPENMP
	#pragma omp barrier
    #endif

    // Move Particles
    #ifdef _OPENMP
	#pragma omp for schedule(dynamic)
    #endif
    for (int i = 0; i < num_parts; ++i) {
        parts[i].move(size);
    }
}



// Command Line Option Processing
int find_arg_idx(int argc, char** argv, const char* option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::atoi(argv[iplace + 1]);
    }

    return default_value;
}

char* find_string_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}

char* find_force_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}

// ==============
// Main Function
// ==============

int main(int argc, char** argv) {
    // Parse Args
    if (find_arg_idx(argc, argv, "-h") >= 0) {
        std::cout << "Options:" << std::endl;
        std::cout << "-h: see this help" << std::endl;
        std::cout << "-n <int>: set number of particles" << std::endl;
        std::cout << "-o <filename>: set the output file name" << std::endl;
        std::cout << "-s <int>: set particle initialization seed" << std::endl;
        std::cout << "-t <int>: set number of threads (working only in parallel mode) [default = 8]" << std::endl;
        return 0;
    }

    // Open Output File
    char* savename = find_string_option(argc, argv, "-o", nullptr);
    if (savename != nullptr) std::cout << "Creating file " << savename << "..." << std::endl;
    std::ofstream fsave(savename);
    if (savename != nullptr) std::cout << "File created." << std::endl;

    //Find force
    char* forcename = find_force_option(argc, argv, "-f", nullptr);
    if (forcename != nullptr) std::cout << "Choosing non default force " <<  forcename << "..." << std::endl;

    const std::unordered_map<std::string, std::shared_ptr<AbstractForce>> fmap =
    {
        {"repulsive", std::make_shared<RepulsiveForce>() },
        {"gravitational", std::make_shared<GravitationalForce>() },
        {"coulomb", std::make_shared<CoulombForce>() },
    };
    
    std::shared_ptr<AbstractForce> force; 
    try {
        force = fmap.at(forcename);      // vector::at throws an out-of-range
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "Default Force chosen "<< '\n';
        force =  fmap.at("repulsive");
    }
    

    // Initialize Particles
    const int num_parts = find_int_arg(argc, argv, "-n", 1000);
    const int part_seed = find_int_arg(argc, argv, "-s", 0);
    const double size = sqrt(density * num_parts);
    const int num_th = find_int_arg(argc, argv, "-t", 8);

    std::vector<particle> parts(num_parts);
    
    std::cout << "Trying to init particles..." << std::endl;
    init_particles(parts, size, part_seed);
    
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    std::cout << "Trying to init simulation..." << std::endl;
    std::cout << "Init simulation ended." << std::endl;

    auto init_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff_1 = init_time - start_time;
    double seconds_1 = diff_1.count();
    std::cout << "initialization Time = " << seconds_1 << " seconds\n";

save(fsave, parts, size);
#ifdef _OPENMP
std::cout << "Available threads: " << std::thread::hardware_concurrency() << "\nRunning "
          << num_th << " thread(s)." <<std::endl;
#pragma omp parallel default(shared) num_threads(num_th)
#endif
    {
        //for nel tempo: non parallelizzare
        for (int step = 0; step < nsteps; ++step) {
            
            simulate_one_step(parts, force, num_parts,size);

            // Save state if necessary
            #ifdef _OPENMP
            #pragma omp master
            #endif
            {
                if (fsave.good() && (step % savefreq) == 0) {
                    save(fsave, parts, size);
                }
                if(step > 0){
                    if (step%10 == 0){
                    fflush(stdout);
                    printf("[ %d% ]\r", (int)(step*100/nsteps));
                    }
                }
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
