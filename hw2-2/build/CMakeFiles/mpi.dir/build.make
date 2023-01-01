# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake

# The command to remove a file.
RM = /u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jellyfish/shared-folder/nbody/hw2-2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jellyfish/shared-folder/nbody/hw2-2/build

# Include any dependencies generated for this target.
include CMakeFiles/mpi.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/mpi.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/mpi.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mpi.dir/flags.make

CMakeFiles/mpi.dir/main.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/main.cpp.o: ../main.cpp
CMakeFiles/mpi.dir/main.cpp.o: CMakeFiles/mpi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mpi.dir/main.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mpi.dir/main.cpp.o -MF CMakeFiles/mpi.dir/main.cpp.o.d -o CMakeFiles/mpi.dir/main.cpp.o -c /home/jellyfish/shared-folder/nbody/hw2-2/main.cpp

CMakeFiles/mpi.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/main.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jellyfish/shared-folder/nbody/hw2-2/main.cpp > CMakeFiles/mpi.dir/main.cpp.i

CMakeFiles/mpi.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/main.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jellyfish/shared-folder/nbody/hw2-2/main.cpp -o CMakeFiles/mpi.dir/main.cpp.s

CMakeFiles/mpi.dir/Particle.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/Particle.cpp.o: ../Particle.cpp
CMakeFiles/mpi.dir/Particle.cpp.o: CMakeFiles/mpi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mpi.dir/Particle.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mpi.dir/Particle.cpp.o -MF CMakeFiles/mpi.dir/Particle.cpp.o.d -o CMakeFiles/mpi.dir/Particle.cpp.o -c /home/jellyfish/shared-folder/nbody/hw2-2/Particle.cpp

CMakeFiles/mpi.dir/Particle.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/Particle.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jellyfish/shared-folder/nbody/hw2-2/Particle.cpp > CMakeFiles/mpi.dir/Particle.cpp.i

CMakeFiles/mpi.dir/Particle.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/Particle.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jellyfish/shared-folder/nbody/hw2-2/Particle.cpp -o CMakeFiles/mpi.dir/Particle.cpp.s

CMakeFiles/mpi.dir/PhysicalForce.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/PhysicalForce.cpp.o: ../PhysicalForce.cpp
CMakeFiles/mpi.dir/PhysicalForce.cpp.o: CMakeFiles/mpi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mpi.dir/PhysicalForce.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mpi.dir/PhysicalForce.cpp.o -MF CMakeFiles/mpi.dir/PhysicalForce.cpp.o.d -o CMakeFiles/mpi.dir/PhysicalForce.cpp.o -c /home/jellyfish/shared-folder/nbody/hw2-2/PhysicalForce.cpp

CMakeFiles/mpi.dir/PhysicalForce.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/PhysicalForce.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jellyfish/shared-folder/nbody/hw2-2/PhysicalForce.cpp > CMakeFiles/mpi.dir/PhysicalForce.cpp.i

CMakeFiles/mpi.dir/PhysicalForce.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/PhysicalForce.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jellyfish/shared-folder/nbody/hw2-2/PhysicalForce.cpp -o CMakeFiles/mpi.dir/PhysicalForce.cpp.s

CMakeFiles/mpi.dir/Find_Arg.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/Find_Arg.cpp.o: ../Find_Arg.cpp
CMakeFiles/mpi.dir/Find_Arg.cpp.o: CMakeFiles/mpi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mpi.dir/Find_Arg.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mpi.dir/Find_Arg.cpp.o -MF CMakeFiles/mpi.dir/Find_Arg.cpp.o.d -o CMakeFiles/mpi.dir/Find_Arg.cpp.o -c /home/jellyfish/shared-folder/nbody/hw2-2/Find_Arg.cpp

CMakeFiles/mpi.dir/Find_Arg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/Find_Arg.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jellyfish/shared-folder/nbody/hw2-2/Find_Arg.cpp > CMakeFiles/mpi.dir/Find_Arg.cpp.i

CMakeFiles/mpi.dir/Find_Arg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/Find_Arg.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jellyfish/shared-folder/nbody/hw2-2/Find_Arg.cpp -o CMakeFiles/mpi.dir/Find_Arg.cpp.s

CMakeFiles/mpi.dir/Output.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/Output.cpp.o: ../Output.cpp
CMakeFiles/mpi.dir/Output.cpp.o: CMakeFiles/mpi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mpi.dir/Output.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mpi.dir/Output.cpp.o -MF CMakeFiles/mpi.dir/Output.cpp.o.d -o CMakeFiles/mpi.dir/Output.cpp.o -c /home/jellyfish/shared-folder/nbody/hw2-2/Output.cpp

CMakeFiles/mpi.dir/Output.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/Output.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jellyfish/shared-folder/nbody/hw2-2/Output.cpp > CMakeFiles/mpi.dir/Output.cpp.i

CMakeFiles/mpi.dir/Output.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/Output.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jellyfish/shared-folder/nbody/hw2-2/Output.cpp -o CMakeFiles/mpi.dir/Output.cpp.s

CMakeFiles/mpi.dir/Simulation.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/Simulation.cpp.o: ../Simulation.cpp
CMakeFiles/mpi.dir/Simulation.cpp.o: CMakeFiles/mpi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mpi.dir/Simulation.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mpi.dir/Simulation.cpp.o -MF CMakeFiles/mpi.dir/Simulation.cpp.o.d -o CMakeFiles/mpi.dir/Simulation.cpp.o -c /home/jellyfish/shared-folder/nbody/hw2-2/Simulation.cpp

CMakeFiles/mpi.dir/Simulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/Simulation.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jellyfish/shared-folder/nbody/hw2-2/Simulation.cpp > CMakeFiles/mpi.dir/Simulation.cpp.i

CMakeFiles/mpi.dir/Simulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/Simulation.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jellyfish/shared-folder/nbody/hw2-2/Simulation.cpp -o CMakeFiles/mpi.dir/Simulation.cpp.s

# Object files for target mpi
mpi_OBJECTS = \
"CMakeFiles/mpi.dir/main.cpp.o" \
"CMakeFiles/mpi.dir/Particle.cpp.o" \
"CMakeFiles/mpi.dir/PhysicalForce.cpp.o" \
"CMakeFiles/mpi.dir/Find_Arg.cpp.o" \
"CMakeFiles/mpi.dir/Output.cpp.o" \
"CMakeFiles/mpi.dir/Simulation.cpp.o"

# External object files for target mpi
mpi_EXTERNAL_OBJECTS =

mpi: CMakeFiles/mpi.dir/main.cpp.o
mpi: CMakeFiles/mpi.dir/Particle.cpp.o
mpi: CMakeFiles/mpi.dir/PhysicalForce.cpp.o
mpi: CMakeFiles/mpi.dir/Find_Arg.cpp.o
mpi: CMakeFiles/mpi.dir/Output.cpp.o
mpi: CMakeFiles/mpi.dir/Simulation.cpp.o
mpi: CMakeFiles/mpi.dir/build.make
mpi: /u/sw/toolchains/gcc-glibc/11.2.0/base/lib/libmpi.so
mpi: CMakeFiles/mpi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable mpi"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mpi.dir/build: mpi
.PHONY : CMakeFiles/mpi.dir/build

CMakeFiles/mpi.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mpi.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mpi.dir/clean

CMakeFiles/mpi.dir/depend:
	cd /home/jellyfish/shared-folder/nbody/hw2-2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jellyfish/shared-folder/nbody/hw2-2 /home/jellyfish/shared-folder/nbody/hw2-2 /home/jellyfish/shared-folder/nbody/hw2-2/build /home/jellyfish/shared-folder/nbody/hw2-2/build /home/jellyfish/shared-folder/nbody/hw2-2/build/CMakeFiles/mpi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mpi.dir/depend

