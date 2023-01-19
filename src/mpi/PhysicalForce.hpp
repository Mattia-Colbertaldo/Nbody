#ifndef HH__PHYSICALFORCE__HH
#define HH__PHYSICALFORCE__HH

#include "Particle.hpp"
#include <vector>

class AbstractForce
{
public:
  AbstractForce(){};
  virtual void force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i, const double mass_n, const double charge_me, const double charge_n) const = 0;
  
  ~AbstractForce(){};
};


class RepulsiveForce : public AbstractForce
{
public:
    
  void force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i,  const double mass_n, const double charge_me, const double charge_n) const override;
  
};
  


class GravitationalForce : public AbstractForce
{
public:
  
  void force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i,  const double mass_n, const double charge_me, const double charge_n) const override;
};


class GravitationalAssistForce : public AbstractForce
{
public:
  
  
  void force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i,  const double mass_n, const double charge_me, const double charge_n) const override;

};

class ProtonForce : public AbstractForce
{
    //equal charged Particles : all are protons
public:
     
  void force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i,  const double mass_n, const double charge_me, const double charge_n) const override;

};



class CoulombForce : public AbstractForce
{
  
public:
   
  void force_application(particle_pos i, particle_pos j, particle_vel_acc& parts_vel_acc_loc_i,  const double mass_n, const double charge_me, const double charge_n) const override;
};




#endif