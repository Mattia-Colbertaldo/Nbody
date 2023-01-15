#include "common.h"
#include "PhysicalForce.hpp"
#include <memory>
    
void RepulsiveForce :: force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i, const double mass_n, const double charge_me, const double charge_n) const {
    // Calculate Distance
    double dx = j.x - i.x;
    double dy = j.y - i.y;
    double dz = j.z - i.z;
    double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);

    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff) return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef = mass_n*(1 - cutoff / r) / r2 ;
    parts_vel_acc_loc_i.ax += coef * dx;
    parts_vel_acc_loc_i.ay += coef * dy;
    parts_vel_acc_loc_i.az += coef * dz;
  
};
  
void GravitationalForce :: force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i, const double mass_n, const double charge_me, const double charge_n) const {
    
    // Calculate Distance
    double dx = j.x - i.x;
    double dy = j.y - i.y;
    double dz = j.z - i.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff) return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef =  (G * mass_n / r2) ;

    parts_vel_acc_loc_i.ax += coef * dx;
    parts_vel_acc_loc_i.ay += coef * dy;
    parts_vel_acc_loc_i.az += coef * dz;
};


void GravitationalAssistForce:: force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i, const double mass_n, const double charge_me, const double charge_n) const {
    // Calculate Distance
    double dx = j.x - i.x;
    double dy = j.y - i.y;
    double dz = j.z - i.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff) return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef;

    // Very simple short-range repulsive force
    if(r2> cutoff * cutoff){
        coef =  G * mass_n / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        coef = -( G * mass_n / r2 ) * 3 ;
    }

    parts_vel_acc_loc_i.ax += coef * dx;
    parts_vel_acc_loc_i.ay += coef * dy;
    parts_vel_acc_loc_i.az += coef * dz;
};



     
void ProtonForce :: force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i, const double mass_n, const double charge_me, const double charge_n) const {
    // Calculate Distance
    double dx = j.x - i.x;
    double dy = j.y - i.y;
    double dz = j.z - i.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff) return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef =  K * proton_charge * proton_charge / r2  ;
    

    parts_vel_acc_loc_i.ax += coef * dx;
    parts_vel_acc_loc_i.ay += coef * dy;
    parts_vel_acc_loc_i.az += coef * dz;
};


   
void CoulombForce :: force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i, const double mass_n, const double charge_me, const double charge_n) const {
    
    // Calculate Distance
    double dx = j.x - i.x;
    double dy = j.y - i.y;
    double dz = j.z - i.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff) return;

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef = std::pow(scale, 2) * K * charge_me * charge_n / r2  ;
    

    parts_vel_acc_loc_i.ax += coef * dx;
    parts_vel_acc_loc_i.ay += coef * dy;
    parts_vel_acc_loc_i.az += coef * dz;
};




