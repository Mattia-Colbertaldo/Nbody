#include "common.cuh"
#include "PhysicalForce.cuh"
#include <cuda.h>

    
void RepulsiveForce :: force_application() const {
    // Calculate Distance
    // double dx = neighbor.x - p.x;
    // double dy = neighbor.y - p.y;
    // double dz = neighbor.z - p.z;
    // double r2 = dx * dx + dy * dy + dz * dz;

    // // Check if the two Particles should interact
    // if (r2 > cutoff * cutoff)
    //     return;

    // r2 = fmax(r2, min_r * min_r);
    // double r = std:: sqrt(r2);

    // // Very simple short-range repulsive force
    // double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
    // p.ax += coef * dx;
    // p.ay += coef * dy;
    // p.az += coef * dz;
  
};
  
void GravitationalForce :: force_application() const {
    // Calculate Distance
    // double dx = neighbor.x - p.x;
    // double dy = neighbor.y - p.y;
    // double dz = neighbor.z - p.z;
    // double r2 = dx * dx + dy * dy + dz * dz;
    // // Check if the two Particles should interact
    // if (r2 > cutoff * cutoff)
    //     return;

    // r2 = fmax(r2, min_r * min_r);
    // double r = std:: sqrt(r2);

    // // Very simple short-range repulsive force
    // double coef =  (G * neighbor.mass / r2) ;

    // p.ax += coef * dx;
    // p.ay += coef * dy;
    // p.az += coef * dz;
};


void GravitationalAssistForce:: force_application() const {
    // Calculate Distance
    // double dx = neighbor.x - p.x;
    // double dy = neighbor.y - p.y;
    // double dz = neighbor.z - p.z;
    // double r2 = dx * dx + dy * dy + dz * dz;
    // // Check if the two Particles should interact
    // if (r2 > cutoff * cutoff)
    //     return;

    // r2 = fmax(r2, min_r * min_r);
    // double r = std:: sqrt(r2);
    // double coef;

    // // Very simple short-range repulsive force
    // if(r2>0.0001){
    //     coef =  G * neighbor.mass / r2 ;
    // }
    // else
    // //gravity-assist : repulsive force
    // {
    //     coef = -( G * neighbor.mass / r2 ) * 3 ;
    // }

    // p.ax += coef * dx;
    // p.ay += coef * dy;
    // p.az += coef * dz;
};



     
void ProtonForce :: force_application() const {
    // Calculate Distance
    // double dx = neighbor.x - p.x;
    // double dy = neighbor.y - p.y;
    // double dz = neighbor.z - p.z;
    // double r2 = dx * dx + dy * dy + dz * dz;
    // // Check if the two Particles should interact
    // if (r2 > cutoff * cutoff)
    //     return;

    // r2 = fmax(r2, min_r * min_r);
    // double r = std:: sqrt(r2);
    // double coef =  K * proton_charge * proton_charge / r2  ;
    

    // p.ax += coef * dx;
    // p.ay += coef * dy;
    // p.az += coef * dz;
};


   
void CoulombForce :: force_application() const {
    
    // // Calculate Distance
    // double dx = neighbor.x - p.x;
    // double dy = neighbor.y - p.y;
    // double dz = neighbor.z - p.z;
    // double r2 = dx * dx + dy * dy + dz * dz;
    // // Check if the two Particles should interact
    // if (r2 > cutoff * cutoff)
    //     return;

    // r2 = fmax(r2, min_r * min_r);
    // double r = std:: sqrt(r2);
    // double coef = std::pow(scale, 2) * K * p.charge * neighbor.charge / r2  ;
    

    // p.ax += coef * dx;
    // p.ay += coef * dy;
    // p.az += coef * dz;
};




