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
#define cutoff   0.01
#define min_r    (cutoff / 100)
#define dt       0.0005



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
        particle(){};
        
        particle(const double x, const double y, const double z,
                const double vx, const double vy, const double vz,
                const float m) : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz),
                ax(0.), ay(0.), az(0.), mass(m){};
        
        void apply_force(const particle& neighbor, const std::string forcename);
        void move(double size);
        

};


class AbstractForce
{
public:

  AbstractForce() = default;
  //! virtual destructor
  virtual ~AbstractForce() = default;

  virtual void  
  force_application(particle& p,const particle& neighbour) const =0;
protected:
    std::string name;    
};


class RepulsiveForce final : public AbstractForce
{
public:
  RepulsiveForce():name( "RepulsiveForce"){};
  virtual ~RepulsiveForce(){} ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual void force_application(particle& p,const particle& neighbor) const override{
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

protected: 
  std::string name;   
};

class GravitationalForce final : public AbstractForce
{
public:
  GravitationalForce():name( "GravitationalForce"){};
  virtual ~GravitationalForce(){} ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual void force_application(particle& p,const particle& neighbor) const override{
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
    constexpr double G= 0.00000000000667;
    // Very simple short-range repulsive force
    if(r2>0.0001){
        coef =  - G * neighbor.mass / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        coef = ( G * neighbor.mass / r2 ) * 3 ;
    }

    p.ax += coef * dx;
    p.ay += coef * dy;
    p.az += coef * dz;
    }
protected:
  std::string name;
};

class CoulombForce final : public AbstractForce
{
    //equal charged particles : all are protons
public:
  CoulombForce():name( "CoulombForce"){};
  virtual ~CoulombForce(){} ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual void force_application(particle& p,const particle& neighbor) const override{
    constexpr double k= 8.98e9;
    constexpr double proton_charge= 1.6e-19;
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
protected:
  std::string name ; 
};



#endif