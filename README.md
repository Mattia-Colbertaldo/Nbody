# N-body simulation
## Overview
This repository contains a versatile N-body simulation program for simulating interactions between particles under different conditions.

* ***Please look at the pdf for an insightful presentation***

This project aims to analyze a wide range of different low level APIs to get the best results for each parallelization option, especially for **CUDA**.

The simulation supports customization of the number of particles, particle seed for reproducibility, and output options. Additionally, it provides parallelization options using OpenMP, MPI, and CUDA for optimized performance. The README provides comprehensive instructions on compilation, execution, and rendering of simulation results.

## Features
* Particle Initialization: Customize the number of particles and set a seed for reproducible initialization.
* Output Control: Specify output file names to record particle positions at each step.
* **Force**: Choose from Gravitational, Gravitational Assist, Anti-Coulomb, Proton and Repulsive forces
* **Collision**: Choose from no collisions, elastic collisions, and inelastic collisions.
* **Parallelization**: Utilize OpenMP, MPI, or CUDA for parallel execution, enhancing simulation performance.
* Rendering: Visualize simulation results using the provided 3drender.py script.

By default, the program runs with 1000 particles. The number of particles can be changed with the "-n" command line parameter:

    Nbody/build> ./serial -n 10000
    Simulation Time = 195.029 seconds for 10000 particles.

If we rerun the program, the initial positions and velocities of the particles will be randomized because the particle seed is unspecified. By default, the particle seed will be unspecified; this can be changed with the "-s" command line parameter:

    Nbody/build> ./serial -s 150
    Simulation Time = 1.45459 seconds for 1000 particles.

This will set the particle seed to 150 which initializes the particles in a reproducible way. You can print the particle positions to a file specified with the "-o" parameter:

    Nbody/build> ./serial -o serial.parts.out
    Simulation Time = 1.78357 seconds for 1000 particles.

This will create a serial.parts.out file with the particle positions after each step listed. You can use the hw2-rendering tool to convert this into a .gif file of your particles.

* Options

    You can use the "-h" command line parameter to print the help menu summarizing the parameter options:

        Nbody/build> ./serial -h
        Options:
        -h: see this help
        -n <int>: set number of particles
        -o <filename>: set the output file name
        -s <int>: set particle initialization seed
        -c [0-1-2]: set collision type (not available in mpi):
            0: no collisions (default)
            1: elastic collision
            2: unelastic collision
        -t <int>: set number of threads (working only in parallel mode) [default = 8]
        -f <string>: choose the force type:
            gravitational
            assist
            coulomb
            proton
            repulsive

* Compiling:
    * Serial, OpenMP and MPI:

            Nbody> mkdir build
            Nbody> cd build
            Nbody/build> cmake ..
    * CUDA:

        Compiling CUDA requires the NVIDIA GPU Computing Toolkit installed, so in the correct machine run:

            Nbody> mkdir build
            Nbody> cd build
            Nbody\build> cmake ..\src\cuda
        
        This will produce a solution file (cuda.sln). Build that solution in Visual Studio, then run with:

            Nbody\build> .\Debug\cuda.exe [args]

* MPI:
    
    * If you are using Docker and the ptrace permissions can't be fixed, run the following command before launching the jobs:

            export OMPI_MCA_btl_vader_single_copy_mechanism=none

        This will disable CMA.

    * Then run with:
        
            mpirun -n <number of processes> mpi [args]

* Rendering:

    * Synopsis:

            python 3drender.py <input_file.txt>

    * Example:
        
            python .\3drender.py ..\build\out.txt
            
 Quoting of the sources:
 the part of the code in the "omp" directory was written on the guidance of https://sites.google.com/lbl.gov/cs267-spr2022/hw-2-1?authuser=0!!br0ken
 
 * Results:

    * Performance:
    
        ![image](https://github.com/Mattia-Colbertaldo/Nbody/assets/100996597/83f33026-6599-4c32-afb9-ecbf5d1e7f84)
    
    * OMP - MPI - CUDA Calculations and Writing times:
      
      (CUDA writing timing includes transfering data from GPU to CPU. A lot of research was done to optimize this operation through the exploitation of various thechniques.)
    
        ![image](https://github.com/Mattia-Colbertaldo/Nbody/assets/100996597/69b71da5-edf2-45b6-8113-40dba77f45e5)

* Here are just a few combinations:
     
     * Gravitational Force - no collisions:

        https://github.com/Mattia-Colbertaldo/Nbody/assets/100996597/adb45f83-b236-424d-8c3b-c98b75670fb5
     
     * Anti - Coulomb (same charge attraction) - no collisions:
    

        https://github.com/Mattia-Colbertaldo/Nbody/assets/100996597/71168cb7-1034-4df1-8989-e36234ed1daf

     * Repulsive - Inelastic collisions:
    

        https://github.com/Mattia-Colbertaldo/Nbody/assets/100996597/84af061a-badf-4bf0-8790-6b5b9eb407a6

    * Proton force - Elastic collisions:
    

        https://github.com/Mattia-Colbertaldo/Nbody/assets/100996597/95ad8c6b-5fc3-4597-a3a6-5c2afca0adb4




     



