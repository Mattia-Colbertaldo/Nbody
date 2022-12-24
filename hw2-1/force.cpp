#include "common.h"
#include <cmath>
#include <vector>

class AbstractForce
{
public:

  AbstractForce() = default;
  //! virtual destructor
  virtual ~AbstractForce() = default;

  virtual particle  
  force_application(particle& p) const
  {
    return p;
  }
  
protected:
    std::string name;    
};


class RepulsiveForce final : public AbstractForce
{
public:
  RepulsiveForce() = default;
  ~RepulsiveForce(){};
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  particle force_application() const override{
    // Calculate Distance
    double dx = neighbor.x - this->x;
    double dy = neighbor.y - this->y;
    double dz = neighbor.z - this->z;
    double r2 = dx * dx + dy * dy + dz * dz;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = neighbor.mass*(1 - cutoff / r) / r2 ;
    this->ax += coef * dx;
    this->ay += coef * dy;
    this->az += coef * dz;
    }

protected:
  std::string name = "RepulsiveForce" ; 
};

class GravitationalForce final : public AbstractForce
{
public:
  GravitationalForce() = default;
  ~GravitationalForce(){};
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  particle force_application() const override{
    // Calculate Distance
    double dx = neighbor.x - me.x;
    double dy = neighbor.y - me.y;
    double dz = neighbor.z - me.z;
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

    this->ax += coef * dx;
    this->ay += coef * dy;
    this->az += coef * dz;
    }
protected:
  std::string name = "GravitationalForce" ; 
};

class CoulombForce final : public AbstractForce
{
    //equal charged particles : all are protons
public:
  CoulombForce() = default;
  ~CoulombForce(){};
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  particle force_application() const override{
    constexpr double k= 8.98 *10^9;
    constexpr double proton_charge= 1.6 *10^(-19);
    // Calculate Distance
    double dx = neighbor.x - me.x;
    double dy = neighbor.y - me.y;
    double dz = neighbor.z - me.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);
    double coef =  k * proton_charge * proton_charge / r2  ;
    

    this->ax += coef * dx;
    this->ay += coef * dy;
    this->az += coef * dz;
    }
protected:
  std::string name = "CoulombForce" ; 
};

