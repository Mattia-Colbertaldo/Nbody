#ifndef __COMMON_H__
#define __COMMON_H__
#include <vector>
#include <random>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <iostream>



// Program Constants
constexpr unsigned int nsteps = 1000;
constexpr unsigned int savefreq = 1;
constexpr double density = 0.0005;
constexpr double cutoff  = 0.01 * 100;
constexpr double min_r   = 0.01 / 1000;
constexpr double dt      = 0.0005;

constexpr double scale = 1e11;

constexpr double G            = 6.67e-11 * scale;
constexpr double K            = 8.98e9 * scale;
constexpr double proton_charge= 1.6e-19 * scale;





#endif