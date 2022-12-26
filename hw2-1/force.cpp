#include "common.h"

#include <cmath>
#include <vector>




void
RepulsiveForce :: force_application(particle& p,const particle& neighbour) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
  
}





void
GravitationalForce :: force_application(particle& p,const particle& neighbour) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    constexpr double G= 6.67 * 10^(-11);

    // Very simple short-range repulsive force
    double coef =  - G * neighbor.mass* / r2 ;

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
}




void
GravitationalAssistForce :: force_application(particle& p,const particle& neighbour) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    constexpr double G= 6.67 * 10^(-11);

    // Very simple short-range repulsive force
    if(r2>0.0001){
        double coef =  - G * neighbor.mass* / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        double coef = ( G * neighbor.mass* / r2 ) * 3 ;
    }

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
}


void
ProtonForce :: force_application(particle& p,const particle& neighbour) const {
    constexpr double k= 8.98 *10^9;
    constexpr double proton_charge= 1.6 *10^(-19);
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef =  k * proton_charge * proton_charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
}



void
CoulombForce :: force_application(particle& p,const particle& neighbour) const {
    constexpr double k= 8.98e9;
    
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef =  k * p.charge * neighbour.charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
}
