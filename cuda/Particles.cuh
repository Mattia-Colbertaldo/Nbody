#ifndef HH__PARTICLE__HH
#define HH__PARTICLE__HH

#include "common.cuh"
#include "PhysicalForce.cuh"
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/transform.h>
#include <thrust/execution_policy.h>
#include <thrust/random.h>
#include <random>

struct saxpy_functor
{
    __host__ __device__
        double operator()(const double& v, const double& x) const {
            return dt * v + x;
        }
};


// __host__ __device__ void saxpy_fast(thrust::device_vector<double>& vx, thrust::device_vector<double>& x)
// {
//    // x <- dt * vx + x
//    thrust::transform(thrust::device, vx.begin(), vx.end(), x.begin(), x.begin(), saxpy_functor());
// }

// struct distance_functor
// {
//     __host__ __device__
//         double operator()(const double& v, const double& x) const {
//         return dt * v + x;
//     }
// };

struct force_functor
{
    thrust::device_vector<double> positions;

     __host__ __device__ force_functor(thrust::device_vector<double>&& positions) : positions(positions) {}

    __host__ __device__ double operator()(thrust::device_vector<double>&& accelerations_ax) const {
        // Calculate Distance
       //  double dx = y - px;
       //  double dy = y - py;
       //  double dz = y - pz;
       //  double r2 = dx * dx + dy * dy + dz * dz;

       //  double coef = 0.;

       //  // Check if the two Particles should interact
       //  if (r2 > cutoff * cutoff){
       //     coef = 0.;
       //  }

       //  else {
       //     r2 = fmax(r2, min_r * min_r);
       //     double r = std:: sqrt(r2);

       //     // Very simple short-range repulsive force
       //     double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
       //  }

       //  p.ax += coef * dx;
       //  p.ay += coef * dy;
       //  p.az += coef * dz;
    }
};


struct Particles {

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

        const int num_parts;
        const double size;

        struct prg
        {
            double a, b;

            __host__ __device__
            prg(double _a, double _b) : a(_a), b(_b) {};

            __host__ __device__
                double operator()(double n) const
                {
                    thrust::default_random_engine rng;
                    thrust::uniform_real_distribution<double> dist(a, b);

                    return (double)dist(rng);
                }
        };

        Particles(const int num_parts, const double size) : num_parts(num_parts), size(size){

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
        };

        // double rand_01()
        // {
        //     std::random_device rd;
        //     std::default_random_engine eng(rd());
        //     std::uniform_real_distribution<double> distr(0., size);
        //     return distr(eng);
        // }

        // double rand_02()
        // {
        //     std::random_device rd;
        //     std::default_random_engine eng(rd());
        //     std::uniform_real_distribution<double> distr(-1.0, 1.0);
        //     return distr(eng);
        // }

        // double rand_03()
        // {
        //     std::random_device rd;
        //     std::default_random_engine eng(rd());
        //     std::uniform_real_distribution<double> distr(-1.0e-19, 1.0e-19);
        //     return distr(eng);
        // }

        __host__ __device__ void Initialize(){


            thrust::counting_iterator<unsigned int> index_sequence_begin(0);

            thrust::transform(x.begin(), x.end(), x.begin(), prg(0.,size));
            thrust::transform(y.begin(), y.end(), y.begin(), prg(0.,size));
            thrust::transform(z.begin(), z.end(), z.begin(), prg(0.,size));

            thrust::transform(vx.begin(), vx.end(), vx.begin(), prg(-1.0, 1.0));
            thrust::transform(vy.begin(), vy.end(), vy.begin(), prg(-1.0, 1.0));
            thrust::transform(vz.begin(), vz.end(), vz.begin(), prg(-1.0, 1.0));

            thrust::transform(masses.begin(), masses.end(), masses.begin(), prg(-1.0, 1.0));
            thrust::transform(charges.begin(), charges.end(), charges.begin(), prg(-1.0e-19, 1.0e-19));


        }

        __host__ __device__ ~Particles(){};

        AbstractForce* force;

        __host__ __device__ void setForce(AbstractForce* force){
            this->force = force;
        };

        
        
        // Particles(const double x, const double y, const double z,
        //         const double vx, const double vy, const double vz,
        //         const double m, const double c) : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz),
        //         ax(0.), ay(0.), az(0.), masses(m), charges(c){};

       __host__ __device__ void move(){
        //  saxpy_fast(ax, vx);
        //  saxpy_fast(ay, vy);
        //  saxpy_fast(az, vz);
        //  saxpy_fast(vx, x);
        //  saxpy_fast(vy, y);
        //  saxpy_fast(vz, z);
       }

        __host__ __device__ void apply_force(AbstractForce* force){
         thrust::fill(ax.begin(), ax.end(), 0.);
         thrust::fill(ax.begin(), ax.end(), 0.);
         thrust::fill(ax.begin(), ax.end(), 0.);
        //  thrust::tranform(thrust::device, ax.begin(), ax.end(), force_functor(positions));
        //  thrust::tranform(thrust::device, ax.begin(), ax.end(), force_functor(positions));
        //  thrust::tranform(thrust::device, ax.begin(), ax.end(), force_functor(positions));

      }
      
       __host__ __device__ void simulate_one_step(){
         //apply_force();
         printf("GPU: Moving particles...");
         move();
         printf("GPU: Particles moved.");
       }

       
};

#endif