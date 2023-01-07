cmake_minimum_required(VERSION 3.17)
project(cuda_test CUDA)

find_package(CUDA)
set(
        CMAKE_CUDA_STANDARD 14
        )

set(DCUDA_NVCC_FLAGS "-arch=sm_60" CACHE STRING "nvcc flags" FORCE)

add_executable(cuda_test main.cu Find_Arg.cu PhysicalForce.cu common.cuh)

set_target_properties(
        cuda_test
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        LINKER_LANGUAGE CXX)