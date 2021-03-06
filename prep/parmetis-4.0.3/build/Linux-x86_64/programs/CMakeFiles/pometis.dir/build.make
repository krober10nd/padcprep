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
include programs/CMakeFiles/pometis.dir/depend.make

# Include the progress variables for this target.
include programs/CMakeFiles/pometis.dir/progress.make

# Include the compile flags for this target's objects.
include programs/CMakeFiles/pometis.dir/flags.make

programs/CMakeFiles/pometis.dir/pometis.c.o: programs/CMakeFiles/pometis.dir/flags.make
programs/CMakeFiles/pometis.dir/pometis.c.o: ../../programs/pometis.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/pometis.dir/pometis.c.o"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/pometis.dir/pometis.c.o   -c /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/pometis.c

programs/CMakeFiles/pometis.dir/pometis.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pometis.dir/pometis.c.i"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -E /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/pometis.c > CMakeFiles/pometis.dir/pometis.c.i

programs/CMakeFiles/pometis.dir/pometis.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pometis.dir/pometis.c.s"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -S /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/pometis.c -o CMakeFiles/pometis.dir/pometis.c.s

programs/CMakeFiles/pometis.dir/pometis.c.o.requires:
.PHONY : programs/CMakeFiles/pometis.dir/pometis.c.o.requires

programs/CMakeFiles/pometis.dir/pometis.c.o.provides: programs/CMakeFiles/pometis.dir/pometis.c.o.requires
	$(MAKE) -f programs/CMakeFiles/pometis.dir/build.make programs/CMakeFiles/pometis.dir/pometis.c.o.provides.build
.PHONY : programs/CMakeFiles/pometis.dir/pometis.c.o.provides

programs/CMakeFiles/pometis.dir/pometis.c.o.provides.build: programs/CMakeFiles/pometis.dir/pometis.c.o

programs/CMakeFiles/pometis.dir/io.c.o: programs/CMakeFiles/pometis.dir/flags.make
programs/CMakeFiles/pometis.dir/io.c.o: ../../programs/io.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/pometis.dir/io.c.o"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/pometis.dir/io.c.o   -c /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/io.c

programs/CMakeFiles/pometis.dir/io.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pometis.dir/io.c.i"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -E /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/io.c > CMakeFiles/pometis.dir/io.c.i

programs/CMakeFiles/pometis.dir/io.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pometis.dir/io.c.s"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && /opt/crc/m/mvapich2/2.1/intel/15.0/mlx/bin/mpicc  $(C_DEFINES) $(C_FLAGS) -S /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs/io.c -o CMakeFiles/pometis.dir/io.c.s

programs/CMakeFiles/pometis.dir/io.c.o.requires:
.PHONY : programs/CMakeFiles/pometis.dir/io.c.o.requires

programs/CMakeFiles/pometis.dir/io.c.o.provides: programs/CMakeFiles/pometis.dir/io.c.o.requires
	$(MAKE) -f programs/CMakeFiles/pometis.dir/build.make programs/CMakeFiles/pometis.dir/io.c.o.provides.build
.PHONY : programs/CMakeFiles/pometis.dir/io.c.o.provides

programs/CMakeFiles/pometis.dir/io.c.o.provides.build: programs/CMakeFiles/pometis.dir/io.c.o

# Object files for target pometis
pometis_OBJECTS = \
"CMakeFiles/pometis.dir/pometis.c.o" \
"CMakeFiles/pometis.dir/io.c.o"

# External object files for target pometis
pometis_EXTERNAL_OBJECTS =

programs/pometis: programs/CMakeFiles/pometis.dir/pometis.c.o
programs/pometis: programs/CMakeFiles/pometis.dir/io.c.o
programs/pometis: programs/CMakeFiles/pometis.dir/build.make
programs/pometis: libparmetis/libparmetis.a
programs/pometis: libmetis/libmetis.a
programs/pometis: programs/CMakeFiles/pometis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable pometis"
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pometis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
programs/CMakeFiles/pometis.dir/build: programs/pometis
.PHONY : programs/CMakeFiles/pometis.dir/build

programs/CMakeFiles/pometis.dir/requires: programs/CMakeFiles/pometis.dir/pometis.c.o.requires
programs/CMakeFiles/pometis.dir/requires: programs/CMakeFiles/pometis.dir/io.c.o.requires
.PHONY : programs/CMakeFiles/pometis.dir/requires

programs/CMakeFiles/pometis.dir/clean:
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs && $(CMAKE_COMMAND) -P CMakeFiles/pometis.dir/cmake_clean.cmake
.PHONY : programs/CMakeFiles/pometis.dir/clean

programs/CMakeFiles/pometis.dir/depend:
	cd /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3 /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/programs /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64 /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs /home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/programs/CMakeFiles/pometis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : programs/CMakeFiles/pometis.dir/depend

