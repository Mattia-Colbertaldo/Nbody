#ifndef HH__OUTPUT__HH
#define HH__OUTPUT__HH
#include <fstream>
#include <vector>
#include "Particles.cuh"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

// =================
// Helper Functions
// =================

/*
    classe che gestisce l’output su file. Ha attributi che identificano il nome del file di output 
    e ogni quanti step fare output su file e const reference ai vettori con le particelle.
    La classe gestisce dietro le quinte l’output, in particolare avrà un metodo save che controlla implicitamente a quale step siamo,
    se a questo step va effettuato l’output e se la risposta è affermativa solo il rank 0 scrive su file. 
*/


class Output
{
public:
    Output(){};
    // I/O routines
    __host__ __device__ void save(const char* fsave, const int num_parts, const Particles particles, const double size, const int& nsteps);

    __host__ __device__ void save_output(const char* fsave, const int num_parts, const int savefreq, const Particles& particles, const int& step, const int& nsteps, const double & size);
};
#endif