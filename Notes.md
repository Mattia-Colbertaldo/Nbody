# Readme and compiling
Readme is clear, MPI and OMP compile with no issues and the executable ends with no segmentation fault. Shows nice scalability.

Some issues when compiling with CUDA (on google colab), had to make some minor changes. When executing the following error appears `ERROR: too many resources requested for launch `, even with one particle and one thread.

# Code
* Code is well organised, good!
* There are some C-style casting. In C++ prefer C++ style casting. They are safer (even if longer to type...).
* I see that you have programmed a parte that shuffles elements. It is ok, but you could have used `std::shuffle` and saved the trouble.
* Sometimes you output using `printf` to have a better control on the output format than just `std::cout`. Yet, since you are using C++, it is nicer to include `cstdio` and use `std::printf`. The syntax is the same (in fact it is just a porting of the C function). 
* Physical parameters are defined as constexpr, this has better performance but a parameter file would be more flexible
* There are couple of old version of CMakeLists files that are not used, they should be deleted
* There are some job files for a cluster that I do not think are used, they should be deleted
* The repo is very heavy (>700MB) due to build files in the history, had to use git-filter-repo to prune it. Put in the git repo only files that are not the product of the compilation process (headers,
source files, documentation).
* Some math functions do not have the `std::` namespace (thus are the one from `math.h`), it is better if all functions come from the standard library (`cmath`). Most (if not all) C standard heaader files have been ported to the standard library with the rule `nameOfCHeader.h` -> `cnameOfCHeader`. 
* Some variables defined inside functions and method that could be const are not const (e.g. `dx` and `dy` in `force_application`). Not a great deal here, but setting const what is meant to be const is a good safety rule.
* `Output` could take as constructor input a `const std::string &` instead of a copy of a string. Again, not a great deal here, but in C++ it is a good rule using `const &` instead of passing-by-value on potentially big objects. You save a copy.
