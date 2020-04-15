# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/parallels/Git/PME/c++/ceres-solver-1.14.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/parallels/Git/PME/c++/ceres-bin

# Include any dependencies generated for this target.
include examples/CMakeFiles/nist.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/nist.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/nist.dir/flags.make

examples/CMakeFiles/nist.dir/nist.cc.o: examples/CMakeFiles/nist.dir/flags.make
examples/CMakeFiles/nist.dir/nist.cc.o: /home/parallels/Git/PME/c++/ceres-solver-1.14.0/examples/nist.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/parallels/Git/PME/c++/ceres-bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/nist.dir/nist.cc.o"
	cd /home/parallels/Git/PME/c++/ceres-bin/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nist.dir/nist.cc.o -c /home/parallels/Git/PME/c++/ceres-solver-1.14.0/examples/nist.cc

examples/CMakeFiles/nist.dir/nist.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nist.dir/nist.cc.i"
	cd /home/parallels/Git/PME/c++/ceres-bin/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/parallels/Git/PME/c++/ceres-solver-1.14.0/examples/nist.cc > CMakeFiles/nist.dir/nist.cc.i

examples/CMakeFiles/nist.dir/nist.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nist.dir/nist.cc.s"
	cd /home/parallels/Git/PME/c++/ceres-bin/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/parallels/Git/PME/c++/ceres-solver-1.14.0/examples/nist.cc -o CMakeFiles/nist.dir/nist.cc.s

examples/CMakeFiles/nist.dir/nist.cc.o.requires:

.PHONY : examples/CMakeFiles/nist.dir/nist.cc.o.requires

examples/CMakeFiles/nist.dir/nist.cc.o.provides: examples/CMakeFiles/nist.dir/nist.cc.o.requires
	$(MAKE) -f examples/CMakeFiles/nist.dir/build.make examples/CMakeFiles/nist.dir/nist.cc.o.provides.build
.PHONY : examples/CMakeFiles/nist.dir/nist.cc.o.provides

examples/CMakeFiles/nist.dir/nist.cc.o.provides.build: examples/CMakeFiles/nist.dir/nist.cc.o


# Object files for target nist
nist_OBJECTS = \
"CMakeFiles/nist.dir/nist.cc.o"

# External object files for target nist
nist_EXTERNAL_OBJECTS =

bin/nist: examples/CMakeFiles/nist.dir/nist.cc.o
bin/nist: examples/CMakeFiles/nist.dir/build.make
bin/nist: lib/libceres.a
bin/nist: /usr/lib/x86_64-linux-gnu/libgflags.so.2.2.1
bin/nist: /usr/lib/x86_64-linux-gnu/libglog.so
bin/nist: /usr/lib/x86_64-linux-gnu/libspqr.so
bin/nist: /usr/lib/x86_64-linux-gnu/libcholmod.so
bin/nist: /usr/lib/x86_64-linux-gnu/libccolamd.so
bin/nist: /usr/lib/x86_64-linux-gnu/libcamd.so
bin/nist: /usr/lib/x86_64-linux-gnu/libcolamd.so
bin/nist: /usr/lib/x86_64-linux-gnu/libamd.so
bin/nist: /usr/lib/x86_64-linux-gnu/liblapack.so
bin/nist: /usr/lib/x86_64-linux-gnu/libf77blas.so
bin/nist: /usr/lib/x86_64-linux-gnu/libatlas.so
bin/nist: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
bin/nist: /usr/lib/x86_64-linux-gnu/librt.so
bin/nist: /usr/lib/x86_64-linux-gnu/libcxsparse.so
bin/nist: /usr/lib/x86_64-linux-gnu/liblapack.so
bin/nist: /usr/lib/x86_64-linux-gnu/libf77blas.so
bin/nist: /usr/lib/x86_64-linux-gnu/libatlas.so
bin/nist: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
bin/nist: /usr/lib/x86_64-linux-gnu/librt.so
bin/nist: /usr/lib/x86_64-linux-gnu/libcxsparse.so
bin/nist: examples/CMakeFiles/nist.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/parallels/Git/PME/c++/ceres-bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/nist"
	cd /home/parallels/Git/PME/c++/ceres-bin/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nist.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/nist.dir/build: bin/nist

.PHONY : examples/CMakeFiles/nist.dir/build

examples/CMakeFiles/nist.dir/requires: examples/CMakeFiles/nist.dir/nist.cc.o.requires

.PHONY : examples/CMakeFiles/nist.dir/requires

examples/CMakeFiles/nist.dir/clean:
	cd /home/parallels/Git/PME/c++/ceres-bin/examples && $(CMAKE_COMMAND) -P CMakeFiles/nist.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/nist.dir/clean

examples/CMakeFiles/nist.dir/depend:
	cd /home/parallels/Git/PME/c++/ceres-bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/parallels/Git/PME/c++/ceres-solver-1.14.0 /home/parallels/Git/PME/c++/ceres-solver-1.14.0/examples /home/parallels/Git/PME/c++/ceres-bin /home/parallels/Git/PME/c++/ceres-bin/examples /home/parallels/Git/PME/c++/ceres-bin/examples/CMakeFiles/nist.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/nist.dir/depend

