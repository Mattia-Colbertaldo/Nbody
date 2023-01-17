cmake_minimum_required(VERSION 3.17)
project(openmp-mpi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# CUDA requires the Toolkit installed, and since we are using Docker to build the other executables, you have to cmake cuda
# in the same machine you have the Toolkit installed: run on the right machine (e.g., from VS) "cmake ../src/cuda"

# find_package(CUDA)

# set( CMAKE_CUDA_STANDARD 14 )

# set(DCUDA_NVCC_FLAGS "-arch=sm_60" "-lineinfo" "-fno-elide-constructors" CACHE STRING "nvcc flags" FORCE)

# We require the GNU compiler for this assignment.
if (NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(Prg-Intel "PrgEnv-intel")
    set(Prg-Clang "PrgEnv-cray")
    # message(WARNING
    #         "Must use GNU Compiler for submission. Make sure you ran:\n"
    #         "module swap ${Prg-${CMAKE_C_COMPILER_ID}} PrgEnv-gnu")
endif ()

add_executable(serial src/openmp/main.cpp src/openmp/Particle.cpp src/openmp/PhysicalForce.cpp src/openmp/Find_Arg.cpp src/openmp/Output.cpp src/openmp/Simulation.cpp src/openmp/common.h)

add_executable(openmp src/openmp/main.cpp src/openmp/Particle.cpp src/openmp/PhysicalForce.cpp src/openmp/Find_Arg.cpp src/openmp/Output.cpp src/openmp/Simulation.cpp src/openmp/common.h)
target_link_libraries(openmp PUBLIC OpenMP::OpenMP_CXX)

add_executable(mpi src/mpi/main.cpp src/mpi/Particle.cpp src/mpi/PhysicalForce.cpp src/mpi/Find_Arg.cpp src/mpi/Output.cpp src/mpi/Simulation.cpp)
target_link_libraries(mpi PRIVATE MPI::MPI_CXX)

# add_executable(cuda src/cuda/main.cu src/cuda/Find_Arg.cu src/cuda/PhysicalForce.cu src/cuda/common.cuh src/cuda/Simulation.cu src/cuda/AllParticles.cu src/cuda/Output.cu)

# set_target_properties(
#         cuda
#         PROPERTIES
#         CUDA_SEPARABLE_COMPILATION ON
#         LINKER_LANGUAGE CXX)

# Copy the job scripts
# configure_file(job-mpi job-mpi COPYONLY)
