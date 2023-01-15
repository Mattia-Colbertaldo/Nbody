#include "Particle.hpp"
#include "Simulation.hpp"
#include <memory>
#include <stdexcept>
#include <cmath>
#include <random>
#include <iostream>


// Particle Initialization
void Simulation::init_particles( const double size, const int part_seed) {
    /*
        input :  1. parts     : vettore di particelle
                 2. masses    : vettore delle masse
                 3. size      : dimensione della particella
                 4. part_seed : seme randomico
    */

    int num_parts = this->parts.size();

    std::random_device rd;
    std::mt19937 gen(part_seed ? part_seed : rd());
    std::mt19937 gen2(part_seed ? part_seed : rd());

    int sx = (int)ceil(std::sqrt((double)num_parts));
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
        std::uniform_real_distribution<double> rand_real(-1.0, 1.0);
        const double vx = rand_real(gen);
        const double vy = rand_real(gen);
        const double vz = rand_real(gen);


        //
        std::uniform_real_distribution<double> rand_mass(0.001, 0.1);
        const double m = rand_mass(gen);

        std::uniform_real_distribution<double> rand_charge(-1.0, 1.0);
        const double charge=rand_charge(gen2) * 1e-19;


        this->parts[i]=Particle(x, y, z, vx, vy, vz, m, charge);
    }

};


    



void Simulation::simulate_one_step( const std::unique_ptr<AbstractForce>& force,const int num_parts,const double size) {
    // Compute Forces
    //int num_parts = parts.size();
    #ifdef _OPENMP
	#pragma omp for schedule(dynamic)
    #endif
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = parts[i].az = 0.;
        for (int j = 0; j < num_parts; ++j) {
            force->force_application(parts[i], parts[j], collision);
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
};
