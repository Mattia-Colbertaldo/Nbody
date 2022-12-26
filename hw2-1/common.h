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


class AbstractForce;

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

  AbstractForce(std::string name) : name(name){};
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
  using AbstractForce::AbstractForce;
  RepulsiveForce():AbstractForce("RepulsiveForce"){};
  virtual ~RepulsiveForce() override = default ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual
  void force_application(particle& p,const particle& neighbor) const override;
protected: 
  std::string name;   
};

class GravitationalForce final : public AbstractForce
{
public:
  using AbstractForce::AbstractForce;
  GravitationalForce():AbstractForce("GravitationalForce"){};
  virtual ~GravitationalForce() override = default ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual
  void force_application(particle& p,const particle& neighbor) const override;
protected:
  std::string name;
};


class GravitationalAssistForce final : public AbstractForce
{
public:
  using AbstractForce::AbstractForce;
  GravitationalAssistForce():AbstractForce("GravitationalAssistForce"){};
 
  virtual ~GravitationalAssistForce() override = default ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual
  void force_application(particle& p,const particle& neighbor) const override;

protected:
  std::string name;
};

class ProtonForce final : public AbstractForce
{
    //equal charged particles : all are protons
public:
  using AbstractForce::AbstractForce;
  ProtonForce():AbstractForce("ProtonForce"){};
 
  virtual ~ProtonForce() override = default ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
  virtual
  void force_application(particle& p,const particle& neighbor) const override;

protected:
  std::string name ; 
};



class CoulombForce final : public AbstractForce
{
    //equal charged particles : all are protons
public:
  using AbstractForce::AbstractForce;
  CoulombForce():AbstractForce("CoulombForce"){};
  virtual ~CoulombForce() override = default ;
  /*!
    The area is positive if vertices are given in
    counterclockwise order
  */
 virtual
  void force_application(particle& p,const particle& neighbor) const override;

protected:
  std::string name ; 
};




#endif