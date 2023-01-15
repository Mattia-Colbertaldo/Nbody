#include "Output.cuh"
#include <iostream>
#include <thrust/device_vector.h>


// =================
// Helper Functions
// =================

/*
    classe che gestisce l’output su file. Ha attributi che identificano il nome del file di output 
    e ogni quanti step fare output su file e const reference ai vettori con le particelle.
    La classe gestisce dietro le quinte l’output, in particolare avrà un metodo save che controlla implicitamente a quale step siamo,
    se a questo step va effettuato l’output e se la risposta è affermativa solo il rank 0 scrive su file. 
*/

void Output::save(std::ofstream& fsave, const std::unique_ptr<AllParticles> & parts, const double size, const int& nsteps){
    
  static bool first = true;

  if (first) {
      fsave << parts->num_parts << " " << size << " " << nsteps << "\n";
      first = false;
  }

  // Opzione 1:
  thrust::copy(parts->x.begin(), parts->x.end(), parts->x_h.begin());
  thrust::copy(parts->y.begin(), parts->y.end(), parts->y_h.begin());
  thrust::copy(parts->z.begin(), parts->z.end(), parts->z_h.begin());
  // Opzione 2:
  // x_h = x;
  // y_h = y;
  // z_h = z;




  cudaDeviceSynchronize();
  for(size_t i = 0; i < parts->num_parts; i++){
    // TODO X_H
        fsave <<  parts->x_h[i] << " " << parts->y_h[i] << " " << parts->z_h[i] << std::endl;
  }
  cudaDeviceSynchronize();
};

void Output::save_output(std::ofstream& fsave, const int savefreq, const std::unique_ptr<AllParticles> & parts , const int& step,  const int& nsteps, const double & size){
    // TODO FIX
    // thrust::copy(x.begin(), x.end(), x_h.begin());
    // thrust::copy(y.begin(), y.end(), y_h.begin());
    // thrust::copy(z.begin(), z.end(), z_h.begin());

    // save(fsave, parts, size, nsteps);
    
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
        std::cout << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;
        std::cout << "Saving: ";
        t = clock();
        for(size_t i = 0; i < parts->num_parts*nsteps; i++){
            if (i%10 == 0){
            fflush(stdout);
            printf("\r");
            printf("Saving: [ %d ]\r", (int)(i*100/(parts->num_parts*nsteps)));
            }
            fsave <<  host_bufferx[i] << " " << host_buffery[i] << " " << host_bufferz[i] << std::endl;
        }
        std::cout << "Saving: " << ((clock() - t)*MS_PER_SEC)/CLOCKS_PER_SEC << " ms" << std::endl;

    }
};
