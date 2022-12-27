#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>
#include <mpi.h>
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




#define G             6.67e-11 * 1000000000000
#define K             8.98e9 * 1000000000000
#define proton_charge 1.6e-19 * 10000000

// Particle Data Structure: used in OPENMP
class particle_pos {
    public:

    particle_pos(){};
    particle_pos(const double x, const double y, const double z) : x(x), y(y), z(z){};
    
    double x;  // Position X
    double y;  // Position Y
    double z;
    
};


class particle_vel_acc {
    public:

    double vx; // Velocity X
    double vy; // Velocity Y
    double vz;
    double ax; // Acceleration X
    double ay; // Acceleration Y
    double az;
    particle_vel_acc(){}; 
    particle_vel_acc(const double vx, const double vy, const double vz) : vx(vx), vy(vy), ax(0.), ay(0.), az(0.){};
    
    void move(particle_pos & me, double size);
};


class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application(particle_vel_acc& particle_vel_acc_loc, particle_pos& p,const particle_pos& neighbor, const float& mass_neigh, const float& charge_me,const float& charge_n  ) const = 0;
  
  virtual ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
    
  void force_application(particle_vel_acc& particle_vel_acc_loc, particle_pos& p,const particle_pos& neighbor, const float& mass_neigh, const float& charge_me,const float& charge_n ) const {
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
    double coef = mass_neigh* (1 - cutoff / r) / r2;
    particle_vel_acc_loc.ax += coef * dx;
    particle_vel_acc_loc.ay += coef * dy;
    particle_vel_acc_loc.az += coef * dz; 
  
};
  
};

class GravitationalForce : public AbstractForce
{
public:
  
  void force_application(particle_vel_acc& particle_vel_acc_loc, particle_pos& p,const particle_pos& neighbor, const float& mass_neigh, const float& charge_me,const float& charge_n ) const {
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
    double coef =  (G * mass_neigh / r2) ;


    particle_vel_acc_loc.ax += coef * dx;
    particle_vel_acc_loc.ay += coef * dy;
    particle_vel_acc_loc.az += coef * dz; 

  };
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application(particle_vel_acc& particle_vel_acc_loc, particle_pos& p,const particle_pos& neighbor, const float& mass_neigh, const float& charge_me,const float& charge_n ) const {

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
        coef =  - G * mass_neigh / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        coef = ( G * mass_neigh / r2 ) * 3 ;
    }


    particle_vel_acc_loc.ax += coef * dx;
    particle_vel_acc_loc.ay += coef * dy;
    particle_vel_acc_loc.az += coef * dz; 
};

};

class ProtonForce : public AbstractForce
{
    //equal charged particles : all are protons
public:
     
  void force_application(particle_vel_acc& particle_vel_acc_loc, particle_pos& p,const particle_pos& neighbor, const float& mass_neigh, const float& charge_me,const float& charge_n ) const {

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


    particle_vel_acc_loc.ax += coef * dx;
    particle_vel_acc_loc.ay += coef * dy;
    particle_vel_acc_loc.az += coef * dz; 

};

};



class CoulombForce : public AbstractForce
{
    //equal charged particles : all are protons
public:
   
  void force_application(particle_vel_acc& particle_vel_acc_loc, particle_pos& p,const particle_pos& neighbor, const float& mass_neigh, const float& charge_me,const float& charge_n ) const {
    double dx = neighbor.x - p.x;
    double dy = neighbor.y - p.y;
    double dz = neighbor.z - p.z;
    double r2 = dx * dx + dy * dy + dz * dz;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;
    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    
    double coef =  K * charge_n * charge_me / r2  ;


    particle_vel_acc_loc.ax += coef * dx;
    particle_vel_acc_loc.ay += coef * dy;
    particle_vel_acc_loc.az += coef * dz; 

    
};

};




#endif