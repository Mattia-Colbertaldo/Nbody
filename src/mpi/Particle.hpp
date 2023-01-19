#ifndef HH__PARTICLE__HH
#define HH__PARTICLE__HH

// Particle Data Structure
struct particle_vel_acc {
        double vx; // Velocity 
        double vy; 
        double vz;
        double ax; // Acceleration
        double ay; 
        double az;
};

struct particle_pos {
        public:
        double x;  // Position 
        double y;  
        double z;

        void move(double size, particle_vel_acc & parts_vel_acc );
};



#endif