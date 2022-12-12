# Nbody

By default, the program runs with 1000 particles. The number of particles can be changed with the "-n" command line parameter:

    student@nid02346:~/hw2-1/build> ./serial -n 10000
    Simulation Time = 195.029 seconds for 10000 particles.

If we rerun the program, the initial positions and velocities of the particles will be randomized because the particle seed is unspecified. By default, the particle seed will be unspecified; this can be changed with the "-s" command line parameter:

    student@nid02346:~/hw2-1/build> ./serial -s 150
    Simulation Time = 1.45459 seconds for 1000 particles.

This will set the particle seed to 150 which initializes the particles in a reproducible way. You can print the particle positions to a file specified with the "-o" parameter:

    student@nid02346:~/hw2-1/build> ./serial -o serial.parts.out
    Simulation Time = 1.78357 seconds for 1000 particles.

This will create a serial.parts.out file with the particle positions after each step listed. You can use the hw2-rendering tool to convert this into a .gif file of your particles.

* Options

    You can use the "-h" command line parameter to print the help menu summarizing the parameter options:

        student@nid02346:~/hw2-1/build> ./serial -h
        Options:
        -h: see this help
        -n <int>: set number of particles
        -o <filename>: set the output file name
        -s <int>: set particle initialization seed
        -t <int>: set number of threads (working only in parallel mode) [default = 8]

* Compiling:

        ../hw2-x> mkdir build
        ../hw2-x> cd build
        ../hw2-x> cmake ..
        ../hw2-x> make

* MPI:
    
    * If you are using docker and the ptrace permissions can't be fixed, run the following command before launching the jobs:

            export OMPI_MCA_btl_vader_single_copy_mechanism=none

        This will disable CMA.

    * Then run with:
        
            mpirun -n <number of processes> mpi [args]

* Checking correctness:

        path/to/correctness/correctness_check.py <file1.txt> <file2.txt>

* Rendering:

    * Synopsis:

            render.py <input_file> <output_file.gif> [cutoff]

    * Example:
        
            ./render.py save.txt image.gif 0.01

    * Requires Python 3 and pillow:

            pip install pillow
