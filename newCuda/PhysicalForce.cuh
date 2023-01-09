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
        
__global__ void  
kernel_no_tiling_force(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts) const override;
  
};
  /*


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


*/

#endif