#include "common.h"
#include "PhysicalForce.hpp"

    
void RepulsiveForce :: force_application(Particle& p,const Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    float m_p = p.mass;
    float m_neigh = neighbor.mass;

    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
           
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
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
  
};
  
void GravitationalForce :: force_application(Particle& p,const Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    float m_p = p.mass;
    float m_neigh = neighbor.mass;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
           
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
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
};


void GravitationalAssistForce:: force_application(Particle& p,const Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    float m_p = p.mass;
    float m_neigh = neighbor.mass;

    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
           
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
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
};



     
void ProtonForce :: force_application(Particle& p,const Particle& neighbor, const int collision) const {
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    float m_p = p.mass;
    float m_neigh = neighbor.mass;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
           
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
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
};


   
void CoulombForce :: force_application(Particle& p,const Particle& neighbor, const int collision) const {
    
    // Calculate Distance
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    float m_p = p.mass;
    float m_neigh = neighbor.mass;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            p.vx = (m_p*p.x + m_neigh*neighbor.x)/(m_p+m_neigh);
            
            p.vy = (m_p*p.y + m_neigh*neighbor.y)/(m_p+m_neigh);
            
            p.vz = (m_p*p.z + m_neigh*neighbor.z)/(m_p+m_neigh);
            
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            p.vx = p.x*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.x*m_neigh/(m_p+m_neigh);
           
            p.vy = p.y*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.y*m_neigh/(m_p+m_neigh);
            
            p.vz = p.z*(m_p-m_neigh)/(m_p + m_neigh) + 2*neighbor.z*m_neigh/(m_p+m_neigh);
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
};




