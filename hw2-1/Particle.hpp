#ifndef HH__PARTICLE__HH
#define HH__PARTICLE__HH

struct Particle {
        Particle(){};
        
        Particle(const double x, const double y, const double z,
                const double vx, const double vy, const double vz,
                const float m, const float c) : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz),
                ax(0.), ay(0.), az(0.), mass(m), charge(c){};
        
        void move(double size);
        
        
        double x;  // Position X
        double y;  // Position Y
        double z;  // Position Z
        double vx; // Velocity X
        double vy; // Velocity Y
        double vz; // Velocity Z
        double ax; // Acceleration X
        double ay; // Acceleration Y
        double az; // Acceleration Z
        float mass;
        float charge;

        

};

#endif