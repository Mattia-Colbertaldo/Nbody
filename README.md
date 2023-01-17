# Nbody

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

* Checking correctness:

        path/to/correctness/check.py <file1.txt> <file2.txt>

* Rendering ( Not working on Docker ) :

    * Synopsis:

            python 3drender.py <input_file.txt>

    * Example:
        
            python .\3drender.py ..\build\out.txt
            
 Quoting of the sources:
 the part of the code in the "omp" directory was written on the guidance of https://sites.google.com/lbl.gov/cs267-spr2022/hw-2-1?authuser=0!!br0ken
