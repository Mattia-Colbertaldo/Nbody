#ifndef HH__PHYSICALFORCE__HH
#define HH__PHYSICALFORCE__HH

#include <cuda.h>

class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application() const = 0;
  
  virtual ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
    
  void force_application() const override;
  
};
  


class GravitationalForce : public AbstractForce
{
public:
  
  void force_application() const override;
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application() const override;

};

class ProtonForce : public AbstractForce
{
    //equal charged Particles& : all are protons
public:
     
  void force_application() const override;

};



class CoulombForce : public AbstractForce
{
    //equal charged Particles& : all are protons
public:
   
  void force_application() const override;
};




#endif