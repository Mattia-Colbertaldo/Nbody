#include "common.h"
#include "PhysicalForce.hpp"
#include <memory>
    
void RepulsiveForce :: force_application(std::vector<particle_pos> parts_pos, std::vector<particle_vel_acc> parts_vel_acc_loc, const double mass_n, const double charge_me, const double charge_n, const int i, const int j, const int collision ) const {
    // Calculate Distance
    double dx = parts_pos[j].x - parts_pos[i].x;
    double dy = parts_pos[j].y - parts_pos[i].y;
    double dz = parts_pos[j].z - parts_pos[i].z;
    double r2 = dx * dx + dy * dy + dz * dz;

    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
            vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

            vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
            vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

            vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
            vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
            vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

            vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
            vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

            vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
            vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef = mass_n*(1 - cutoff / r) / r2 ;
    parts_vel_acc_loc[i].ax += coef * dx;
    parts_vel_acc_loc[i].ay += coef * dy;
    parts_vel_acc_loc[i].az += coef * dz;
  
};
  
void GravitationalForce :: force_application(std::vector<particle_pos> parts_pos, std::vector<particle_vel_acc> parts_vel_acc_loc, const double mass_n, const double charge_me, const double charge_n, const int i, const int j, const int collision ) const {
    // Calculate Distance
    double dx = parts_pos[j].x - parts_pos[i].x;
    double dy = parts_pos[j].y - parts_pos[i].y;
    double dz = parts_pos[j].z - parts_pos[i].z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
            vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

            vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
            vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

            vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
            vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
            vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

            vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
            vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

            vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
            vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);

    // Very simple short-range repulsive force
    double coef =  (G * mass_n / r2) ;

    parts_vel_acc_loc[i].ax += coef * dx;
    parts_vel_acc_loc[i].ay += coef * dy;
    parts_vel_acc_loc[i].az += coef * dz;
};


void GravitationalAssistForce:: force_application(std::vector<particle_pos> parts_pos, std::vector<particle_vel_acc> parts_vel_acc_loc, const double mass_n, const double charge_me, const double charge_n, const int i, const int j, const int collision ) const {
    // Calculate Distance
    double dx = parts_pos[j].x - parts_pos[i].x;
    double dy = parts_pos[j].y - parts_pos[i].y;
    double dz = parts_pos[j].z - parts_pos[i].z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;
    
    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
            vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

            vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
            vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

            vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
            vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
            vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

            vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
            vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

            vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
            vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef;

    // Very simple short-range repulsive force
    if(r2>0.0001){
        coef =  G * mass_n / r2 ;
    }
    else
    //gravity-assist : repulsive force
    {
        coef = -( G * mass_n / r2 ) * 3 ;
    }

    parts_vel_acc_loc[i].ax += coef * dx;
    parts_vel_acc_loc[i].ay += coef * dy;
    parts_vel_acc_loc[i].az += coef * dz;
};



     
void ProtonForce :: force_application(std::vector<particle_pos> parts_pos, std::vector<particle_vel_acc> parts_vel_acc_loc, const double mass_n, const double charge_me, const double charge_n, const int i, const int j, const int collision ) const {
    // Calculate Distance
    double dx = parts_pos[j].x - parts_pos[i].x;
    double dy = parts_pos[j].y - parts_pos[i].y;
    double dz = parts_pos[j].z - parts_pos[i].z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;
    
    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
            vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

            vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
            vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

            vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
            vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
            vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

            vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
            vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

            vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
            vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef =  K * proton_charge * proton_charge / r2  ;
    

    parts_vel_acc_loc[i].ax += coef * dx;
    parts_vel_acc_loc[i].ay += coef * dy;
    parts_vel_acc_loc[i].az += coef * dz;
};


   
void CoulombForce :: force_application(std::vector<particle_pos> parts_pos, std::vector<particle_vel_acc> parts_vel_acc_loc, const double mass_n, const double charge_me, const double charge_n, const int i, const int j, const int collision ) const {
    
    // Calculate Distance
    double dx = parts_pos[j].x - parts_pos[i].x;
    double dy = parts_pos[j].y - parts_pos[i].y;
    double dz = parts_pos[j].z - parts_pos[i].z;
    double r2 = dx * dx + dy * dy + dz * dz;
    // Check if the two Particles should interact
    if (r2 > cutoff * cutoff)
        return;

    if(r2 < min_r*min_r){
        if(collision > 0){
            
            // TODO ARGUMENT
            if(collision== 1){
            // URTO ANELASTICO:
            vx[thx] = (mx*tile1_vx[threadIdx.x] + my*tile2_vx[threadIdx.y])/(mx+my);
            vx[thy] = (my*tile2_vx[threadIdx.y] + mx*tile1_vx[threadIdx.x])/(my+mx);

            vy[thx] = (mx*tile1_vy[threadIdx.x] + my*tile2_vy[threadIdx.y])/(mx+my);
            vy[thy] = (my*tile2_vy[threadIdx.y] + mx*tile1_vy[threadIdx.x])/(my+mx);

            vz[thx] = (mx*tile1_vz[threadIdx.x] + my*tile2_vz[threadIdx.y])/(mx+my);
            vz[thy] = (my*tile2_vz[threadIdx.y] + mx*tile1_vz[threadIdx.x])/(my+mx);
            }
            // "unelastic" collision
            else if(collision== 2){
            // URTO ELASTICO
            vx[thx] = tile1_vx[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vx[threadIdx.y]*my/(mx+my);
            vx[thy] = tile2_vx[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vx[threadIdx.x]*mx/(mx+my);

            vy[thx] = tile1_vy[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vy[threadIdx.y]*my/(mx+my);
            vy[thy] = tile2_vy[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vy[threadIdx.x]*mx/(mx+my);

            vz[thx] = tile1_vz[threadIdx.x]*(mx-my)/(mx + my) + 2*tile2_vz[threadIdx.y]*my/(mx+my);
            vz[thy] = tile2_vz[threadIdx.y]*(my-mx)/(mx + my) + 2*tile1_vz[threadIdx.x]*mx/(mx+my);
            }
        }
        return;
    }

    r2 = fmax(r2, min_r * min_r);
    double r = std:: sqrt(r2);
    double coef = std::pow(scale, 2) * K * charge_me * charge_n / r2  ;
    

    parts_vel_acc_loc[i].ax += coef * dx;
    parts_vel_acc_loc[i].ay += coef * dy;
    parts_vel_acc_loc[i].az += coef * dz;
};




