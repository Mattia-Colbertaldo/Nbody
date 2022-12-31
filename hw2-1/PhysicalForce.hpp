#ifndef HH__PHYSICALFORCE__HH
#define HH__PHYSICALFORCE__HH

#include "Particle.hpp"

class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application(Particle& p,const Particle& neighbor) const = 0;
  
  virtual ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
    
  void force_application(Particle& p,const Particle& neighbor) const override;
  
};
  


class GravitationalForce : public AbstractForce
{
public:
  
  void force_application(Particle& p,const Particle& neighbor) const override;
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application(Particle& p,const Particle& neighbor) const override;

};

class ProtonForce : public AbstractForce
{
    //equal charged Particles : all are protons
public:
     
  void force_application(Particle& p,const Particle& neighbor) const override;

};



class CoulombForce : public AbstractForce
{
    //equal charged Particles : all are protons
public:
   
  void force_application(Particle& p,const Particle& neighbor) const override;
};




#endif