#ifndef HH__PARTICLE__HH
#define HH__PARTICLE__HH

// Particle Data Structure: used in OPENMP
struct particle_pos {
        double x;  // Position X
        double y;  // Position Y
        double z;
};

struct particle_vel_acc {
        double vx; // Velocity X
        double vy; // Velocity Y
        double vz;
        double ax; // Acceleration X
        double ay; // Acceleration Y
        double az;
};


struct Particle {
        public:
        Particle(const double & x, const double & y, const double & z){
           part_pos.x=x;   
           part_pos.y=y;    
           part_pos.z=z;     
        } 
        
        void move(double size);
        
        particle_pos part_pos;
        particle_vel_acc part_vel_acc;

        double mass;
        double charge;

};



#endif