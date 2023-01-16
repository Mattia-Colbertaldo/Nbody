#ifndef HH__FIND_ARG__HH
#define HH__FIND_ARG__HH
#include <cmath>
#include <memory>
#include "PhysicalForce.hpp"

class Find_Arg
{
    public:
        Find_Arg(int argc, char** argv):argc(argc), argv(argv), cl(argc,argv){};
        
        int find_int_arg( const std::string type_of_find, const int default_value) const;
        std::string find_string_arg(const std::string type_of_find, const std::string default_value) const;
        std::unique_ptr<AbstractForce> find_force(const std::string forcename) const;
        
    private:
        int argc;
        char** argv;
        GetPot cl;
};
#endif
