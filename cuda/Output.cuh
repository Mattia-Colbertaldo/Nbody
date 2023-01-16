#ifndef HH__OUTPUT__HH
#define HH__OUTPUT__HH
#include <fstream>
#include <vector>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "common.cuh"
#include "AllParticles.cuh"
#include <sstream>
// =================
// Helper Functions
// =================

/*
    classe che gestisce l’output su file. Ha attributi che identificano il nome del file di output 
    e ogni quanti step fare output su file e const reference ai vettori con le particelle.
    La classe gestisce dietro le quinte l’output, in particolare avrà un metodo save che controlla implicitamente a quale step siamo,
    se a questo step va effettuato l’output e se la risposta è affermativa solo il rank 0 scrive su file. 
*/

using namespace common_h;

class Output
{
public:
    Output(const int num_parts, const std::string filename) : num_parts(num_parts), filename(filename){
        bufferx.resize(num_parts*nsteps);
        buffery.resize(num_parts*nsteps);
        bufferz.resize(num_parts*nsteps);
        
        host_bufferx.resize(num_parts*nsteps);
        host_buffery.resize(num_parts*nsteps);
        host_bufferz.resize(num_parts*nsteps);


    };
    // I/O routines
    void save( const std::unique_ptr<AllParticles> & parts, const double size);

    void save_output( const std::unique_ptr<AllParticles> & parts , const int& step, const double & size);
 

    private:
        int num_parts;
        std::string filename;
        std::ostringstream strstream;
        //CAMBIAMENTO: SONO PRIVATI
        thrust::device_vector<double> bufferx;
        thrust::device_vector<double> buffery;
        thrust::device_vector<double> bufferz;
        thrust::host_vector<double> host_bufferx;
        thrust::host_vector<double> host_buffery;
        thrust::host_vector<double> host_bufferz;

};
#endif