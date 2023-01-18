#include "Find_Arg.hpp"
#include <iostream>
#include <memory>

/*/
template<typename T>
std::unique_ptr<T> make_unique() {
    return std::unique_ptr<T>(new T());
}
*/

void printHelp()
{
    std::cout << "Options:" << std::endl;
    std::cout << "-h: see this help" << std::endl;
    std::cout << "-n <int>: set number of particles" << std::endl;
    std::cout << "-o <filename>: set the output file name" << std::endl;
    std::cout << "-s <int>: set particle initialization seed" << std::endl;
    std::cout << "-t <int>: set number of threads (working only in parallel mode) [default = 8]" << std::endl;
    std::cout << "-f <int>: set force: default, repulsive, gravitational, assist, proton, coulomb" << std::endl;
    return ;
};

std::string Find_Arg::find_string_option ( const std::string type_of_find, const std::string default_value) {
    
    std::string result = cl.follow(default_value.c_str(), type_of_find.c_str());
    return result;
};




int Find_Arg :: find_int_arg( const std::string type_of_find, const int default_value) {
    if(type_of_find== "-h"){
        if(cl.search("-h")) {
            printHelp();
            return 0;
        }
    }
    
};


std::unique_ptr<AbstractForce> Find_Arg::find_force(const std::string forcename) 
{
    std::unique_ptr<AbstractForce> force;   
    if(forcename.compare("gravitational")==0){
        force = std::make_unique<GravitationalForce>();
    }
    
    else if(forcename.compare("assist")==0){
        force = std::make_unique<GravitationalAssistForce>();
    }
    
    else if(forcename.compare("proton")==0){
        
        force = std::make_unique<ProtonForce>();
    }
    
    else if(forcename.compare("coulomb")==0){
        force = std::make_unique<CoulombForce>();
    }

    else {
        force = std::make_unique<RepulsiveForce>();
    }
    return force;
    
};


