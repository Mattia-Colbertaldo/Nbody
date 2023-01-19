#include "Output.hpp"
#include <iostream>


using namespace common_h;
// I/O routines
void Output :: save( const std::vector<Particle>& parts, const double size) {
    int num_parts = parts.size();
    
    static bool first = true;

    if (first) {
        this->strstream << num_parts << " " << size << " " << nsteps << "\n";
        first = false;
    }

    for (int i = 0; i < num_parts; ++i) {
        this->strstream << parts[i].x << " " << parts[i].y << " " << parts[i].z << "\n";
    }

};


void Output :: save_output( const std::vector<Particle>& parts , const int& step, const double & size)
{
    if (strstream.good() && (step % savefreq) == 0)
    {
        save(parts, size);
        if(step == nsteps-1){
            long t = clock();
            std::ofstream(filename) << this->strstream.str();
            std::cout << "Writing: " << ((clock() - t)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;
            return;
        }
    }
    if(step > 0 && step < nsteps - 1){
        if (step%10 == 0){
        fflush(stdout);
        printf("[ %d% ]\r", (int)(step*100/nsteps));
        }
    }
};