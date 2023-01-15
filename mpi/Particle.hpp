#ifndef HH__PARTICLE__HH
#define HH__PARTICLE__HH

// Particle Data Structure
struct particle_vel_acc {
        double vx; // Velocity X
        double vy; // Velocity Y
        double vz;
        double ax; // Acceleration X
        double ay; // Acceleration Y
        double az;
};

struct particle_pos {
        public:
        double x;  // Position X
        double y;  // Position Y
        double z;

        void move(double size, particle_vel_acc & parts_vel_acc );
};



#endif