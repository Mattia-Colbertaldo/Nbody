#ifndef HH__OUTPUT__HH
#define HH__OUTPUT__HH
#include <fstream>
#include <vector>
#include <sstream>
#include "Particle.hpp"
#include "common.h"
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
    Output(std::string filename) : filename(filename){};
    // I/O routines
    void save(const std::vector<Particle>& parts, const double size);

    void save_output(const std::vector<Particle>& parts , const int& step, const double & size);

private:
    std::ostringstream strstream;
    std::string filename;
};
#endif