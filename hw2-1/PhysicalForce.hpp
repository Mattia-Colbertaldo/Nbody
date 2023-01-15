#ifndef HH__PHYSICALFORCE__HH
#define HH__PHYSICALFORCE__HH

#include "Particle.hpp"

class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application(Particle& p,const Particle& neighbor, const int collision) const = 0;
  
  ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
    
  void force_application(Particle& p,const Particle& neighbor, const int collision) const override;
  
};
  


class GravitationalForce : public AbstractForce
{
public:
  
  void force_application(Particle& p,const Particle& neighbor, const int collision) const override;
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application(Particle& p,const Particle& neighbor, const int collision) const override;

};

class ProtonForce : public AbstractForce
{
    //equal charged Particles : all are protons
public:
     
  void force_application(Particle& p,const Particle& neighbor, const int collision) const override;

};



class CoulombForce : public AbstractForce
{
    //equal charged Particles : all are protons
public:
   
  void force_application(Particle& p,const Particle& neighbor, const int collision) const override;
};




#endif