#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>
#include <random>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <iostream>


// Program Constants
#define nsteps   1000
#define savefreq 1
#define density  0.0005
#define cutoff   0.01 * 100
#define min_r    0.01 / 1000
#define dt       0.0005

#define scale 1e11

#define G             6.67e-11 * scale
#define K             8.98e9 * scale
#define proton_charge 1.6e-19 * scale



// Particle Data Structure: used in OPENMP

class particle {
    public:
        
        double x;  // Position X
        double y;  // Position Y
        double z;  // Position Z
        double vx; // Velocity X
        double vy; // Velocity Y
        double vz; // Velocity Z
        double ax; // Acceleration X
        double ay; // Acceleration Y
        double az; // Acceleration Z
        float mass;
        float charge;
        particle(){};
        
        particle(const double x, const double y, const double z,
                const double vx, const double vy, const double vz,
                const float m, const float c) : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz),
                ax(0.), ay(0.), az(0.), mass(m), charge(c){};
        
        void move(double size);
        

};


class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application(particle& p,const particle& neighbor) const = 0;
  
  virtual ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
    
  void force_application(particle& p,const particle& neighbor) const {
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
  
};
  
};

class GravitationalForce : public AbstractForce
{
public:
  
  void force_application(particle& p,const particle& neighbor) const {
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
    double coef =  (G * neighbor.mass / r2) ;

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
  };
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application(particle& p,const particle& neighbor) const {
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

};

class ProtonForce : public AbstractForce
{
    //equal charged particles : all are protons
public:
     
  void force_application(particle& p,const particle& neighbor) const {
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
    double coef =  K * proton_charge * proton_charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
};

};



class CoulombForce : public AbstractForce
{
    //equal charged particles : all are protons
public:
   
  void force_application(particle& p,const particle& neighbor) const {
    
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
    double coef = std::pow(scale, 2) * K * p.charge * neighbor.charge / r2  ;
    

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
};

};




#endif