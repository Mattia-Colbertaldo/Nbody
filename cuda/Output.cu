#include "Output.cuh"
#include <iostream>
#include <cstdio>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/execution_policy.h>


// =================
// Helper Functions
// =================

/*
    classe che gestisce l’output su file. Ha attributi che identificano il nome del file di output 
    e ogni quanti step fare output su file e const reference ai vettori con le particelle.
    La classe gestisce dietro le quinte l’output, in particolare avrà un metodo save che controlla implicitamente a quale step siamo,
    se a questo step va effettuato l’output e se la risposta è affermativa solo il rank 0 scrive su file. 
*/

struct functor
{
  __host__ __device__
  void operator()(double val)
  {
      printf("%f ", val);
  }
};

// I/O routines
__host__ __device__ void Output :: save(const char* charname, const int num_parts, const Particles particles, const double size, const int& nsteps) {
    
    printf("  GPU: Saving...\n");
    printf("  GPU: %d %f %d", num_parts, size, nsteps);
    static bool first = true;

    if (first) {
        printf("%d %f %d", num_parts, size, nsteps);
        first = false;
    }

    for (int i = 0; i < num_parts; ++i) {
        printf("%f %f %f", (double)(particles.x[i]), (double)particles.y[i], (double)particles.z[i]);
        //thrust::for_each(thrust::device, particles.positions.get<0>().begin(),particles.positions.get<0>().end(),functor());
    }

    // fsave << std::endl;
};


__host__ __device__ void Output :: save_output(const char* charname, const int num_parts, const int savefreq, const Particles& particles, const int& step, const int& nsteps, const double & size)
{
    if ((step % savefreq) == 0)
    {
        save(charname, num_parts, particles, size, nsteps);
    }
    // if(step > 0){
    //     if (step%10 == 0){
    //     fflush(stdout);
    //     printf("[ %d% ]\r", (int)(step*100/nsteps));
    //     }
    // }
};
