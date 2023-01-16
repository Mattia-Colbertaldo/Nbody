#ifndef HH__FIND_ARG__HH
#define HH__FIND_ARG__HH
#include <cmath>
#include <memory>
#include "PhysicalForce.hpp"
#include "GetPot"

class Find_Arg
{
    public:
        Find_Arg(int argc, char** argv):argc(argc), argv(argv), cl(argc,argv){
            //GetPot temp(argc, argv);
            //cl=temp;
        };
        
        int find_int_arg( const std::string type_of_find, const int default_value);
        std::string find_string_option( const std::string type_of_find, const std::string default_value);
        std::unique_ptr<AbstractForce> find_force(const std::string forcename);
        
    private:
        int argc;
        char** argv;
        GetPot cl;

};
#endif
