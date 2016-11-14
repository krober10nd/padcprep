# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /afs/crc.nd.edu/x86_64_linux/cmake/3.2.2/bin/cmake

# The command to remove a file.
RM = /afs/crc.nd.edu/x86_64_linux/cmake/3.2.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64

# Include any dependencies generated for this target.
include programs/CMakeFiles/parmetis_prog.dir/depend.make

# Include the progress variables for this target.
include programs/CMakeFiles/parmetis_prog.dir/progress.make

# Include the compile flags for this target's objects.
include programs/CMakeFiles/parmetis_prog.dir/flags.make

programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o: programs/CMakeFiles/parmetis_prog.dir/flags.make
programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o: ../../programs/parmetis.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/parmetis_prog.dir/parmetis.c.o   -c /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/parmetis.c

programs/CMakeFiles/parmetis_prog.dir/parmetis.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/parmetis_prog.dir/parmetis.c.i"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -E /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/parmetis.c > CMakeFiles/parmetis_prog.dir/parmetis.c.i

programs/CMakeFiles/parmetis_prog.dir/parmetis.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/parmetis_prog.dir/parmetis.c.s"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -S /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/parmetis.c -o CMakeFiles/parmetis_prog.dir/parmetis.c.s

programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.requires:
.PHONY : programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.requires

programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.provides: programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.requires
	$(MAKE) -f programs/CMakeFiles/parmetis_prog.dir/build.make programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.provides.build
.PHONY : programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.provides

programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.provides.build: programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o

programs/CMakeFiles/parmetis_prog.dir/io.c.o: programs/CMakeFiles/parmetis_prog.dir/flags.make
programs/CMakeFiles/parmetis_prog.dir/io.c.o: ../../programs/io.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/parmetis_prog.dir/io.c.o"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/parmetis_prog.dir/io.c.o   -c /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/io.c

programs/CMakeFiles/parmetis_prog.dir/io.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/parmetis_prog.dir/io.c.i"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -E /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/io.c > CMakeFiles/parmetis_prog.dir/io.c.i

programs/CMakeFiles/parmetis_prog.dir/io.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/parmetis_prog.dir/io.c.s"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -S /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/io.c -o CMakeFiles/parmetis_prog.dir/io.c.s

programs/CMakeFiles/parmetis_prog.dir/io.c.o.requires:
.PHONY : programs/CMakeFiles/parmetis_prog.dir/io.c.o.requires

programs/CMakeFiles/parmetis_prog.dir/io.c.o.provides: programs/CMakeFiles/parmetis_prog.dir/io.c.o.requires
	$(MAKE) -f programs/CMakeFiles/parmetis_prog.dir/build.make programs/CMakeFiles/parmetis_prog.dir/io.c.o.provides.build
.PHONY : programs/CMakeFiles/parmetis_prog.dir/io.c.o.provides

programs/CMakeFiles/parmetis_prog.dir/io.c.o.provides.build: programs/CMakeFiles/parmetis_prog.dir/io.c.o

programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o: programs/CMakeFiles/parmetis_prog.dir/flags.make
programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o: ../../programs/adaptgraph.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/parmetis_prog.dir/adaptgraph.c.o   -c /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/adaptgraph.c

programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/parmetis_prog.dir/adaptgraph.c.i"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -E /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/adaptgraph.c > CMakeFiles/parmetis_prog.dir/adaptgraph.c.i

programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/parmetis_prog.dir/adaptgraph.c.s"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -S /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/adaptgraph.c -o CMakeFiles/parmetis_prog.dir/adaptgraph.c.s

programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.requires:
.PHONY : programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.requires

programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.provides: programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.requires
	$(MAKE) -f programs/CMakeFiles/parmetis_prog.dir/build.make programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.provides.build
.PHONY : programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.provides

programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.provides.build: programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o

# Object files for target parmetis_prog
parmetis_prog_OBJECTS = \
"CMakeFiles/parmetis_prog.dir/parmetis.c.o" \
"CMakeFiles/parmetis_prog.dir/io.c.o" \
"CMakeFiles/parmetis_prog.dir/adaptgraph.c.o"

# External object files for target parmetis_prog
parmetis_prog_EXTERNAL_OBJECTS =

programs/parmetis: programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o
programs/parmetis: programs/CMakeFiles/parmetis_prog.dir/io.c.o
programs/parmetis: programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o
programs/parmetis: programs/CMakeFiles/parmetis_prog.dir/build.make
programs/parmetis: libparmetis/libparmetis.a
programs/parmetis: libmetis/libmetis.a
programs/parmetis: programs/CMakeFiles/parmetis_prog.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable parmetis"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/parmetis_prog.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
programs/CMakeFiles/parmetis_prog.dir/build: programs/parmetis
.PHONY : programs/CMakeFiles/parmetis_prog.dir/build

programs/CMakeFiles/parmetis_prog.dir/requires: programs/CMakeFiles/parmetis_prog.dir/parmetis.c.o.requires
programs/CMakeFiles/parmetis_prog.dir/requires: programs/CMakeFiles/parmetis_prog.dir/io.c.o.requires
programs/CMakeFiles/parmetis_prog.dir/requires: programs/CMakeFiles/parmetis_prog.dir/adaptgraph.c.o.requires
.PHONY : programs/CMakeFiles/parmetis_prog.dir/requires

programs/CMakeFiles/parmetis_prog.dir/clean:
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && $(CMAKE_COMMAND) -P CMakeFiles/parmetis_prog.dir/cmake_clean.cmake
.PHONY : programs/CMakeFiles/parmetis_prog.dir/clean

programs/CMakeFiles/parmetis_prog.dir/depend:
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3 /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64 /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs/CMakeFiles/parmetis_prog.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : programs/CMakeFiles/parmetis_prog.dir/depend
