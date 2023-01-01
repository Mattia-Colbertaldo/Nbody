#include "Find_Arg.hpp"
#include <iostream>
#include <memory>


// Command Line Option Processing
int find_arg_idx(int argc, char** argv, const std::string option) {
    for (int i = 1; i < argc; ++i) {
        std::string str(argv[i]) ;
        if (option.compare(str) == 0) {
            return i;
        }
    }
    return -1;
};

std::string find_string_option(int argc, char** argv, const std::string option, std::string default_value) {
    int iplace = find_arg_idx(argc, argv, option);
    std::cout << option << " " << iplace  << std::endl;

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
};

std::string find_force_option(int argc, char** argv, const std::string option, std::string default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
};

std::string Find_Arg :: find_string_arg(const std::string type_of_find, const std::string option){
    if("-o"==type_of_find)
    {
        return find_string_option(this->argc, this->argv, type_of_find, option);                
    }
    else if("-f"==type_of_find)
    {
        return find_force_option(this->argc, this->argv, type_of_find, option);                
    }
};

int Find_Arg :: find_int_arg( const std::string type_of_find, const int default_value){
    int iplace = find_arg_idx(this->argc, this->argv, type_of_find);
    if("-h"==type_of_find) return iplace;
    else if(iplace >= 0 && iplace < this->argc - 1) {
        return std::atoi(this->argv[iplace + 1]);
    }

    return default_value;
};


std::shared_ptr<AbstractForce> Find_Arg::find_force(const std::string forcename)
{
    std::shared_ptr<AbstractForce> force;   
    if(forcename.compare("gravitational")==0){
        std::cout << "Gravitational force chosen." << std::endl;
        force = std::make_shared<GravitationalForce>();
    }
    
    else if(forcename.compare("assist")==0){
        std::cout << "Gravitational Assist force chosen." << std::endl;
        force = std::make_shared<GravitationalAssistForce>();
    }
    
    else if(forcename.compare("proton")==0){
        
        std::cout << "Proton force chosen." << std::endl;
        force = std::make_shared<ProtonForce>();
    }
    
    else if(forcename.compare("coulomb")==0){
        std::cout << "Coulomb force chosen." << std::endl;
        force = std::make_shared<CoulombForce>();
    }

    else {
        std::cout << "Repulsive force chosen." << std::endl;
        force = std::make_shared<RepulsiveForce>();
    }
    return force;
    
};


