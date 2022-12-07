#include "common.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <mpi.h>

// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream& fsave, std::vector<particle_mpi>& parts, int num_parts, double size) {
    //int num_parts = parts.size();
    static bool first = true;

    if (first) {
        fsave << num_parts << " " << size << "\n";
        first = false;
    }

    for (int i = 0; i < num_parts; ++i) {
        fsave << parts[i].x << " " << parts[i].y << "\n";
    }

    fsave << std::endl;
}

// Particle Initialization
void init_particles(std::vector<particle_mpi>& parts, std::vector<float>& masses, int num_parts, double size,int part_seed) {
    //int num_parts = parts.size();
    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());

    int sx = (int)ceil(sqrt((double)num_parts));
    int sy = (num_parts + sx - 1) / sx;

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

        // Distribute particles evenly to ensure proper spacing
        parts[i].x = size * (1. + (k % sx)) / (1 + sx);
        parts[i].y = size * (1. + (k / sx)) / (1 + sy);

        /*
        /*
        // Assign random velocities within a bound
        std::uniform_real_distribution<float> rand_real(-1.0, 1.0);
        parts[i].vx = rand_real(gen);
        parts[i].vy = rand_real(gen);
        */

        // Assing random mass
        std::uniform_real_distribution<float> rand_mass(0.001, 0.1);
        float m = rand_mass(gen);
        masses.emplace_back(m);
    }
    //std::cout << masses << std::endl;
    //std::cout << masses << std::endl;
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

    // Initialize Particles
    int num_parts = find_int_arg(argc, argv, "-n", 1000);
    int part_seed = find_int_arg(argc, argv, "-s", 0);
    double size = sqrt(density * num_parts);
    int num_th = find_int_arg(argc, argv, "-t", 8);

   

    std::vector<particle_mpi> parts(num_parts);
    std::vector<float> masses(num_parts);
    std::cout << "Trying to init particles..." << std::endl;
    init_particles(parts, masses, num_parts, size, part_seed);
    int rank, mpi_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
    // Algorithm
    auto start_time = std::chrono::steady_clock::now();

    std::cout << "Trying to init simulation..." << std::endl;
    init_simulation(parts, masses, num_parts, size);
    std::cout << "Init simulation ended." << std::endl;

    auto init_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff_1 = init_time - start_time;
    double seconds_1 = diff_1.count();
    std::cout << "initialization Time = " << seconds_1 << " seconds\n";




    
    //for nel tempo: non parallelizzare
    for (int step = 0; step < nsteps; ++step) {
        std::cout << "0\n";
        simulate_one_step(parts, masses, num_parts, size);
        std::cout << "One step simulated.\n";

        // Save state if necessary
        if(rank==0)
        {/*
            if (fsave.good() && (step % savefreq) == 0) {
                gather_for_save(parts, masses, num_parts, rank, size );
            }
            */
            if(step > 0){
                if (step%10 == 0){
                fflush(stdout);
                printf("[ %d% ]\r", (int)(step*100/nsteps));
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
    //delete[] parts;
    MPI_Finalize();
}
