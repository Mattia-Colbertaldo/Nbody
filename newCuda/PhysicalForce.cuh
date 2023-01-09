#ifndef HH__PHYSICALFORCE__HH
#define HH__PHYSICALFORCE__HH

#include <cuda.h>

class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, const dim3 grid_sizes, const dim3 block_sizes) const = 0;
  
  virtual ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
        
void  
force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, const dim3 grid_sizes, const dim3 block_sizes) const override;
  
};
  


class GravitationalForce : public AbstractForce
{
public:
  
  void force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, const dim3 grid_sizes, const dim3 block_sizes) const override;
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, const dim3 grid_sizes, const dim3 block_sizes) const override;

};

class ProtonForce : public AbstractForce
{
    //equal charged Particles& : all are protons
public:
     
  void force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, const dim3 grid_sizes, const dim3 block_sizes) const override;

};



class CoulombForce : public AbstractForce
{
    //equal charged Particles& : all are protons
public:
   
  void force_application(double* x, double* y, double* z, double* vx, double* vy, double* vz,
                        double* ax, double* ay, double* az, const double* masses, const double* charges, const int num_parts, const dim3 grid_sizes, const dim3 block_sizes) const override;
};




#endif