#ifndef HH__ALLPARTICLES__HH
#define HH__ALLPARTICLES__HH

struct AllParticles {
        AllParticles(){};
        
        AllParticles(const int num_parts, AbstractForce* force) : num_parts(num_parts), force(force){
            this->size = std::sqrt(density * num_parts);
            
            thrust::default_random_engine rng;

            //Option 1
            pos.reserve(num_parts);
            vel.reserve(num_parts);
            acc.reserve(num_parts);
            // initialize_0_size(pos, rng, size);
            // initialize_11(vel, rng);
            thrust::fill(acc.begin(), acc.end(), make_double3(0.0, 0.0, 0.0));
            dpos = thrust::raw_pointer_cast(pos.data());
            dvel = thrust::raw_pointer_cast(vel.data());
            dacc = thrust::raw_pointer_cast(acc.data());
            // Option 2
            x_h.reserve(num_parts);
            y_h.reserve(num_parts);
            z_h.reserve(num_parts);
            x.reserve(num_parts);
            y.reserve(num_parts);
            z.reserve(num_parts);
            vx.reserve(num_parts);
            vy.reserve(num_parts);
            vz.reserve(num_parts);
            ax.reserve(num_parts);
            ay.reserve(num_parts);
            az.reserve(num_parts);
            masses.reserve(num_parts);
            charges.reserve(num_parts);

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
        };
        
        void move(const double size, double* dx, double* dy, double* dz,
                    double* dvx, double* dvy, double* dvz,
                    double* dax, double* day, double* daz, const double size, const int num_parts , 
                    dim3 grid_sizes, const dim3 block_sizes);
        
        
        public:
        // Option 1
        thrust::device_vector<double3> pos;
        thrust::device_vector<double3> vel;
        thrust::device_vector<double3> acc;
        double3* dpos;
        double3* dvel;
        double3* dacc;
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
        const int num_parts;
        double size;
        AbstractForce* force;

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