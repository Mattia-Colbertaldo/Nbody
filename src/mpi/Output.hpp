#ifndef HH__OUTPUT__HH
#define HH__OUTPUT__HH
#include <fstream>
#include <vector>
#include <sstream>
#include "common.h"
#include "Simulation.hpp"

class Output
{
public:
    Output(std::string filename) : filename(filename){};
    // I/O routines
    void save(const std::vector<particle_pos>& parts, const double size);

    void save_output(const std::vector<particle_pos>& parts , const int& step, const double & size);

private:
    std::ostringstream strstream;
    std::string filename;
};
#endif


