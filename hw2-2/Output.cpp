#include "Output.hpp"
#include <iostream>


// =================
// Helper Functions
// =================

/*
    classe che gestisce l’output su file. Ha attributi che identificano il nome del file di output 
    e ogni quanti step fare output su file e const reference ai vettori con le particelle.
    La classe gestisce dietro le quinte l’output, in particolare avrà un metodo save che controlla implicitamente a quale step siamo,
    se a questo step va effettuato l’output e se la risposta è affermativa solo il rank 0 scrive su file. 
*/

// I/O routines
void Output :: save(std::ofstream& fsave, const std::vector<particle_pos>& parts, const double size, const int& nsteps) {
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
};


void Output :: save_output(std::ofstream& fsave, const int savefreq, const std::vector<particle_pos>& parts , const int& step,  const int& nsteps, const double & size)
{
    if (fsave.good() && (step % savefreq) == 0)
    {
        save(fsave, parts, size, nsteps);
    }
    if(step > 0){
        if (step%10 == 0){
        fflush(stdout);
        printf("[ %d% ]\r", (int)(step*100/nsteps));
        }
    }
};
