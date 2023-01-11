#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <thrust/device_vector.h>

/*
 * --- Using Thrust vectors in kernel functions ---
 * This example demonstrates how we can use device vectors (vectors that live
 * in GPU) as arguments to a kernel launched by the CPU.
 *
 * Every time you push a value in the thrust::device_vector it is expensive as
 * it needs to send it individually to GPU. So in this case it is a not very
 * good way to initialize a particle system, but it is good for an example. It
 * would be better to initialize all its (vector) members all at once.
 *
 * The device_vectors are members of the particle system. Remember that Thrust
 * is simply an abstraction that works on CPU. So we cannot send them to a
 * kernel. However what we can do, is to get the GPU pointer to the first
 * element of the vector and send it to the kernel. This is absolutely great.
 *
 * Compilation instructions (in NCCA labs):
 *     nvcc -g -I/usr/include/cuda -x cu main.cpp 
 * 
 * NOTE: Thrust lives under /usr/include/cuda/thrust so with one -I we can find
 * both by prefixing <thrust/...> or just <...> (for cuda headers).
 */

// This is the kernel that is launched from CPU and GPU runs it for each cell
__global__ void kernel(float *positions, float *velocities)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
    positions[index] += velocities[index];
    // positions[index] = index; // use this one for debugging the index
}

// The particle system stores the point positions and velocities in
// thrust::vector (not to be confused with a cartesian space vector). Cool thing
// of vectors is that they store the elements contiguously in memory. If we had
// 3 particles in the system the positions thrust::device_vector would look like
//     positions: [particle0posx, particle0posy, particle0posz, particle1posx...]
// Same applies for the velocities. That is a sensible way of storing data if we
// want later OpenGL to render them. I will send an example on how to share a
// Vertex attribute buffer with CUDA later.
class ParticleSystem
{
  public:
    ParticleSystem() : num(0) {};
    void birth_particle();       // pushes back one more particle data to the device_vectors
    void advance_particles();    // launches the kernel that adds velocity to positions
    thrust::device_vector<float> positions;   
    thrust::device_vector<float> velocities;
    unsigned int num;            // number of particles
};

int main(void)
{
  int numParticles = 1024;
  int numSteps = 1;    // understand it as "frames", how many steps in time

  ParticleSystem ps;
  
  // Send new positions and velocities to GPU
  for (int i = 0; i < numParticles; i++)
    ps.birth_particle();
  
  // Particle data lives in GPU now so we call the kernel on them few times!
  // This is great! As we don't have to be retrieving and re-sending, Thrust
  // functionality shines in this step. Great framework.
  for (int i = 0; i < numSteps; i++)
    ps.advance_particles();

  std::cout << std::fixed;
	std::cout << std::setprecision(2);

  // This is gonna be an expensive way of accessing the positions, as for each
  // call to the ::operator[]() we are fetching the item from GPU to CPU. The 
  // right way would be to copy the device_vector into a host_vector like this:
  //    thrust::host_vector<float> host_vector = ps.positions
  // That would do it paralelly, it would transfer all the items from the device
  // vector into the host_vector in a parallel way, but to keep it simply in the
  // code I will not be using it.
  for (int i = 0; i < ps.num; i++)
    std::cout << ps.positions[3*i] << " " << ps.positions[3*i+1] << " " << ps.positions[3*i+2] << std::endl;

  return 0;
}

void ParticleSystem::birth_particle()
{
  positions.push_back(2.0f);
  positions.push_back(2.5f);
  positions.push_back(2.0f);

  velocities.push_back(2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f);
  velocities.push_back(2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f);
  velocities.push_back(2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f);

  num += 1;
}

void ParticleSystem::advance_particles()
{
  // As we cannot send device vectors to the kernel (as device_vector is at
  // the end of the day a GPU structure abstraction in CPU) we have to get the
  // pointer in GPU memory in order for the kernel to know where to start 
  // reading the float arrays from.
  float* d_positions =  thrust::raw_pointer_cast(&positions[0]);
  float* d_velocities = thrust::raw_pointer_cast(&velocities[0]);

  /* This is the way I structured my blocks and threads. I fixed the amount of
   * threads per block to 1024. So to get the amount of blocks we just get the
   * total number of elements in positions and divide it by 1024. We add one in
   * case the division leaves remainder.
   *
   * ┌──────────────────────grid─┬of─blocks─────────────────┬──────────
   * │     block_of_threads      │     block_of_threads     │  
   * │ ┌───┬───┬───────┬──────┐  │ ┌───┬───┬───────┬──────┐ │
   * │ │ 0 │ 1 │ [...] │ 1023 │  │ │ 0 │ 1 │ [...] │ 1023 │ │   ...
   * │ └───┴───┴───────┴──────┘  │ └───┴───┴───────┴──────┘ │
   * └───────────────────────────┴──────────────────────────┴──────────
   */
  unsigned int num_of_elements_per_array = 3 * num;
  unsigned int block_size = 1024;
  unsigned int grid_size = num_of_elements_per_array / block_size + 1;

  std::cout << "Num of elements per array: " << num_of_elements_per_array << std::endl;
  std::cout << "Num of blocks in grid: " << grid_size << std::endl;
  std::cout << "Num of threads per block: " << block_size << std::endl;

  // Launch the kernel! As you can see we are not copying memory from CPU to GPU
  // as you would normally do with cudaMemcpy(), as we don't need to! The
  // vectors live in GPU already so we just need to know where they start (GPU
  // pointer) and pass it to the kernel. No need to copy back, we can read from
  // the device vector with the ::operator[]() i.e. positions[2] and that would
  // do all the memory copying for us!
  kernel<<<grid_size,block_size>>>(d_positions, d_velocities);
}