#ifndef HH_FORCE_HH
#define HH_FORCE_HH
#include "common.h"
#include "particle.hpp"
#include <memory>
#include <cmath>
#include <vector>
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
  virtual ~RepulsiveForce(){} override;
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
  virtual ~GravitationalForce(){} override;
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
protected:
  std::string name;
};

class CoulombForce final : public AbstractForce
{
    //equal charged particles : all are protons
public:
  CoulombForce():name( "CoulombForce"){};
  virtual ~CoulombForce(){} override;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual void force_application(particle& p,const particle& neighbor) const override{
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
protected:
  std::string name ; 
};

#endif