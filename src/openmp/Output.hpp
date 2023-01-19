#ifndef HH__OUTPUT__HH
#define HH__OUTPUT__HH
#include <fstream>
#include <vector>
#include <sstream>
#include "Particle.hpp"
#include "common.h"


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