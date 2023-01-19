#include "common.h"
#include "PhysicalForce.hpp"

using namespace common_h;

void RepulsiveForce :: force_application(Particle& p, Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
    float m_p = p.mass;
    float m_neigh = neighbor.mass;

    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            if(collision== 1){
            // inelastic collision:
            #ifdef _OPENMP
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            #endif
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            neighbor.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            neighbor.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            neighbor.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            }
            
            else if(collision== 2){
            // elasic collision
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);

            neighbor.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            neighbor.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            neighbor.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
            
            }
            
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    neighbor.az -= coef * dz;
  
};
  
void GravitationalForce :: force_application(Particle& p, Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
    float m_p = p.mass;
    float m_neigh = neighbor.mass;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            if(collision== 1){
            // inelastic collision:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            neighbor.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            neighbor.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            neighbor.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            }
            
            else if(collision== 2){
            // elasic collision
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);

            neighbor.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            neighbor.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            neighbor.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
            
            }
        
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef =  (G * neighbor.mass / r2) ;

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    neighbor.az -= coef * dz;
};


void GravitationalAssistForce:: force_application(Particle& p, Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
    float m_p = p.mass;
    float m_neigh = neighbor.mass;

    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            if(collision== 1){
            // inelastic collision:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            neighbor.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            neighbor.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            neighbor.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            }
            
            else if(collision== 2){
            // elasic collision
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);

            neighbor.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            neighbor.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            neighbor.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
            
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef;

    // Very simple short-range repulsive force
    if(r2>0.0001){
        coef =  G * neighbor.mass / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        coef = -( G * neighbor.mass / r2 ) * 3 ;
    }

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    neighbor.az -= coef * dz;
};



     
void ProtonForce :: force_application(Particle& p, Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
    float m_p = p.mass;
    float m_neigh = neighbor.mass;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            if(collision== 1){
            // inelastic collision:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            neighbor.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            neighbor.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            neighbor.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            }
            
            else if(collision== 2){
            // elasic collision
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);

            neighbor.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            neighbor.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            neighbor.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
            
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef =  K * proton_charge * proton_charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    neighbor.az -= coef * dz;
};


   
void CoulombForce :: force_application(Particle& p, Particle& neighbor, const int collision) const {
    
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
    float m_p = p.mass;
    float m_neigh = neighbor.mass;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            if(collision== 1){
            // inelastic collision:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            neighbor.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            neighbor.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            neighbor.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            }
            
            else if(collision== 2){
            // elasic collision
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);

            neighbor.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
            neighbor.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            neighbor.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
            
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef = std::pow(scale, 2) * K * p.charge * neighbor.charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    neighbor.az -= coef * dz;
};




