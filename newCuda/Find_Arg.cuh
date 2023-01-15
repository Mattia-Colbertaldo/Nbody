#ifndef HH__FIND_ARG__HH
#define HH__FIND_ARG__HH
#include <cmath>
#include <memory>
#include "PhysicalForce.cuh"
#include <string>

class Find_Arg
{
    public:
        Find_Arg(int argc, char** argv):argc(argc), argv(argv){};
        
        int find_int_arg( const std::string option, const int default_value) const;
        std::string find_string_arg(const std::string type_of_find, const std::string option) const;
        std::shared_ptr<AbstractForce> find_force(const std::string forcename) const;
        
    private:
        int argc;
        char** argv;
};
#endif
