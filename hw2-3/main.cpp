#include "common.h"

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

#define MATRIX_SIZE 1000*1000
#define CPU_VECTOR_SIZE 1024
#define MASK_WIDTH 221
#define TILE_WIDTH 221




// =================
// Helper Functions
// =================

// I/O routines
void save(std::ofstream& fsave, const std::vector<particle>& parts, const double size) {
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
void init_particles(std::vector<particle>& parts, const double size, const int part_seed) {
    /*
        input :  1. parts     : vettore di particelle
                 2. masses    : vettore delle masse
                 3. size      : dimensione della particella
                 4. part_seed : seme randomico
    */

    int num_parts = parts.size();

    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());
    std::mt19937 gen2(part_seed ? part_seed : rd());

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

        std::uniform_real_distribution<float> rand_charge(-1.0, 1.0);
        const double charge=rand_charge(gen2) * 1e-19;


        //
        std::uniform_real_distribution<float> rand_mass(0.001, 0.1);
        const float m = rand_mass(gen);
        parts[i]=particle(x, y, z, vx, vy, vz, m, charge);
    }

}

AbstractForce* Find_force(const char* forcename)
{
    AbstractForce* force;
    if(strcmp(forcename, "gravitational")==0){
        GravitationalForce* f = new GravitationalForce();
        std::cout << "Gravitational force chosen." << std::endl;
        force = f;
    }
    
    else if(strcmp(forcename, "assist")==0){
        GravitationalAssistForce* f = new GravitationalAssistForce();
        std::cout << "Gravitational Assist force chosen." << std::endl;
        force = f;
    }
    
    else if(strcmp(forcename, "proton")==0){
        ProtonForce* f = new ProtonForce();
        std::cout << "Proton force chosen." << std::endl;
        force = f;
    }
    
    else if(strcmp(forcename, "coulomb")==0){
        CoulombForce* f = new CoulombForce();
        std::cout << "Coulomb force chosen." << std::endl;
        force = f;
    }

    else {
        RepulsiveForce* f = new RepulsiveForce();
        std::cout << "Repulsive force chosen." << std::endl;
        force = f;
    }
    return force;
    
}


void simulate_one_step(std::vector<particle>& parts, const AbstractForce& force,const int num_parts,const double size) {
    
    int thx=threadIdx.x;
    
    int pos = blockIdx.x * TILE_WIDTH + thx;
    
    if (pos<num_parts) {
      // Compute Forces
      //int num_parts = parts.size();
      
      parts[pos].ax = parts[pos].ay = parts[pos].az = 0.;
      for (int j = 0; j < num_parts; ++j) {
          force.force_application(parts[pos], parts[j]);
      }
      
    }

    parts[pos].move(size);
    
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
        std::cout << "-f <int>: set force: default, repulsive, gravitational, assist, proton, coulomb" << std::endl;
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
    else{
        std::string def="default";
        forcename= &def[0];
        std::cout << "Choosing default force..." << std::endl;;
    }

    AbstractForce* force= Find_force(forcename);
    std::vector<particle> parts(num_parts);

  
    
    

    // Initialize Particles
    const int num_parts = find_int_arg(argc, argv, "-n", 1000);
    const int part_seed = find_int_arg(argc, argv, "-s", 0);
    const double size = sqrt(density * num_parts);
    const int num_th = find_int_arg(argc, argv, "-t", 8);


    unsigned char *parts;
    parts = cudaMallocManaged((void **) &parts, sizeof(particle)*(num_parts));
    force = cudaMallocManaged((void **) &force, sizeof(AbstractForce));
    num_parts = cudaMallocManaged((void **) &num_parts, sizeof(int)));
    size = cudaMallocManaged((void **) &size, sizeof(double));

    std::cout << "Trying to init particles..." << std::endl;
    init_particles(parts, size, part_seed);

    int deviceId;
    cudaGetDevice(&deviceId);
    cudaMemPrefetchAsync(parts,sizeof(particle)*(num_parts),deviceId);
    cudaMemPrefetchAsync(force,sizeof(AbstractForce),deviceId);
    cudaMemPrefetchAsync(num_parts,sizeof(int),deviceId);
    cudaMemPrefetchAsync(size,sizeof(double),deviceId);
      
    // some events to count the execution time
    //clock_t st, end;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    std::cout << "Trying to init simulation..." << std::endl;
    std::cout << "Init simulation ended." << std::endl;



    cudaEventRecord(start, 0);

    //for nel tempo: non parallelizzare
    for (int step = 0; step < nsteps; ++step) {
        
        simulate_one_stepl<<<num_parts, TILE_WIDTH>>>(parts, *force, num_parts, size );

        // Save state if necessary
        
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
        
    
    cudaThreadSynchronize();

    // time counting terminate
    
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    // compute time elapsed on GPU computing
    cudaEventElapsedTime(&naive_gpu_elapsed_time_ms, start, stop);
    

    save(fsave, parts, size);

    // Finalize
    std::cout << "Simulation Time = " << naive_gpu_elapsed_time_ms << " seconds for " << num_parts <<
     " particles and " << nsteps << " steps.\n";
    fsave.close();

    // free memory
    cudaFree(parts);
    cudaFree(force);
    cudaFree(num_parts);
    cudaFree(size);

}
