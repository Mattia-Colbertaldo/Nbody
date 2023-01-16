#ifndef HH__ALLPARTICLES__HH
#define HH__ALLPARTICLES__HH



#include "common.cuh"
#include "PhysicalForce.cuh"

#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <sm_60_atomic_functions.h>
#include <fstream>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/random.h>
#include <thrust/transform.h>
#include <thrust/execution_policy.h>

#include <thrust/zip_function.h>
#include <thrust/iterator/zip_iterator.h>

struct AllParticles {
    
        public:

        AllParticles(const int num_parts, std::shared_ptr<AbstractForce>  force) : num_parts(num_parts){

            this->force = std :: move(force); 
            this->size = std::sqrt(density * num_parts);
            
            thrust::default_random_engine rng;

            // Option 2
            x_h.resize(num_parts);
            y_h.resize(num_parts);
            z_h.resize(num_parts);
            x.resize(num_parts);
            y.resize(num_parts);
            z.resize(num_parts);
            vx.resize(num_parts);
            vy.resize(num_parts);
            vz.resize(num_parts);
            ax.resize(num_parts);
            ay.resize(num_parts);
            az.resize(num_parts);
            masses.resize(num_parts);
            charges.resize(num_parts);

            // Copying to host vectors (Optional)
            
            // Block and Grid dimensions
            th_per_block = 32;
            block_sizes.x = 32;
            block_sizes.y = 32;
            // block_sizes.x = 1024;
            // block_sizes.y = 1024;
            // block_sizes.z = 0;
            // grid_sizes(ceil(num_parts/th_per_block), ceil(num_parts/th_per_block));
            grid_sizes.x = 1;
            grid_sizes.y = 1;
            // grid_sizes.x = 2147483647/1024;
            // grid_sizes.y = 65535/1024;
            // grid_sizes.z = 0;

            dx = thrust::raw_pointer_cast(x.data());
            dy = thrust::raw_pointer_cast(y.data());
            dz = thrust::raw_pointer_cast(z.data());
            hx = thrust::raw_pointer_cast(x.data());
            hy = thrust::raw_pointer_cast(y.data());
            hz = thrust::raw_pointer_cast(z.data());
            dvx = thrust::raw_pointer_cast(vx.data());
            dvy = thrust::raw_pointer_cast(vy.data());
            dvz = thrust::raw_pointer_cast(vz.data());
            dax = thrust::raw_pointer_cast(ax.data());
            day = thrust::raw_pointer_cast(ay.data());
            daz = thrust::raw_pointer_cast(az.data());
            dmasses = thrust::raw_pointer_cast(masses.data());
            dcharges = thrust::raw_pointer_cast(charges.data());
        };
        
        void init();
        void move();
        void ResetAccelerations();
        
        // Option 2
        thrust::host_vector<double> x_h;
        thrust::host_vector<double> y_h;
        thrust::host_vector<double> z_h;
        thrust::device_vector<double> x;
        thrust::device_vector<double> y;
        thrust::device_vector<double> z;
        thrust::device_vector<double> vx;
        thrust::device_vector<double> vy;
        thrust::device_vector<double> vz;
        thrust::device_vector<double> ax;
        thrust::device_vector<double> ay;
        thrust::device_vector<double> az;
        thrust::device_vector<double> masses;
        thrust::device_vector<double> charges;
        double* dx;
        double* dy;
        double* dz;
        double* hx;
        double* hy;
        double* hz;
        double* dvx;
        double* dvy;
        double* dvz;
        double* dax;
        double* day;
        double* daz;
        double* dmasses;
        double* dcharges;

        private:
        //CAMBIAMENTO 2: private
            const int num_parts;
            double size;
            std::shared_ptr<AbstractForce> force;

            int th_per_block = 10;
            dim3 block_sizes;
            // block_sizes.x;
            // block_sizes.y;
            // block_sizes.z;
            // dim3 grid_sizes;
            dim3 grid_sizes;
            // grid_sizes.x;
            // grid_sizes.y;
            // grid_sizes.z;

        

};




#endif