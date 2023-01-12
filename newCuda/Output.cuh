#ifndef HH__OUTPUT__HH
#define HH__OUTPUT__HH
#include <fstream>
#include <vector>
#include "AllParticles.cuh"
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
    void save(std::ofstream& fsave, const std::shared_ptr<AllParticles> & parts, const double size, const int& nsteps);

    void save_output(std::ofstream& fsave, const int savefreq, const std::shared_ptr<AllParticles> & parts , const int& step,  const int& nsteps, const double & size);
};
#endif