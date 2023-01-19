#include "Output.cuh"
#include <iostream>
#include <thrust/device_vector.h>
#include <sstream>


using namespace common_h;
void Output::save( const std::unique_ptr<AllParticles> & parts, const double size){
    
  static bool first = true;

  if (first) {
      std::cout << "first!" << std::endl;
      strstream << num_parts << " " << size << " " << nsteps << "\n";
      first = false;
  }

  thrust::copy(parts->x.begin(), parts->x.end(), parts->x_h.begin());
  thrust::copy(parts->y.begin(), parts->y.end(), parts->y_h.begin());
  thrust::copy(parts->z.begin(), parts->z.end(), parts->z_h.begin());

  cudaDeviceSynchronize();
  for(size_t i = 0; i < num_parts; i++){
    strstream <<  parts->x_h[i] << " " << parts->y_h[i] << " " << parts->z_h[i] << "\n";
  }
  cudaDeviceSynchronize();
};

void Output::save_output(const std::unique_ptr<AllParticles> & parts , const int& step, const double & size){
    
    thrust::copy(parts->x.begin(), parts->x.end(), &bufferx[num_parts*step]);
    thrust::copy(parts->y.begin(), parts->y.end(), &buffery[num_parts*step]);
    thrust::copy(parts->z.begin(), parts->z.end(), &bufferz[num_parts*step]);

    cudaDeviceSynchronize();
    if(step > 0){
        if (step%10 == 0){
        fflush(stdout);
        printf("[ %d ]\r", (int)(step*100/nsteps));
        }
    }
    if(step == nsteps - 1){
        std::cout << "Retrieving data from GPU: ";
        long t = clock();
        thrust::copy(bufferx.begin(), bufferx.end(), host_bufferx.begin());
        thrust::copy(buffery.begin(), buffery.end(), host_buffery.begin());
        thrust::copy(bufferz.begin(), bufferz.end(), host_bufferz.begin());
        std::cout << ((clock() - t)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;
        std::cout << "Saving: ";
        t = clock();
        for(size_t i = 0; i < num_parts*nsteps; i++){
            if (i%10 == 0){
            fflush(stdout);
            printf("\r");
            printf("Saving: [ %d ]\r", (int)(i*100/(num_parts*nsteps)));
            }
            strstream <<  host_bufferx[i] << " " << host_buffery[i] << " " << host_bufferz[i] << "\n";
        }
        std::ofstream(filename) << strstream.str();
        std::cout << "Saving: " << ((clock() - t)*1000)/CLOCKS_PER_SEC << " ms" << std::endl;

    }
};
