#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>
#include <random>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <cuda.h>

namespace common_h{

    constexpr float MS_PER_SEC = 1000.0f;
    // constexpr int grid_size = num_parts / block_size + 1;

    constexpr unsigned int BLOCK_DIM = 32;

    // Program Constants
    constexpr unsigned int nsteps = 1000;
    constexpr unsigned int savefreq = 1;
    constexpr double density = 0.0005;
    constexpr double cutoff  = 0.01 * 100;
    constexpr double min_r   = 0.01 / 10;
    constexpr double dt      = 0.0005;
    constexpr double max_velocity = 1;
    constexpr int common_block_size = 32;

    constexpr double scale = 1e10;

    constexpr double G            = 6.67e-11 * scale;
    constexpr double K            = 8.98e9 * scale;
    constexpr double proton_charge= 1.6e-19 * scale;
}


#endif