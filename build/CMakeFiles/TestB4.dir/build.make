# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

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

# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_COMMAND = /home/alexsta1993/miniforge3/envs/alexandros/bin/cmake

# The command to remove a file.
RM = /home/alexsta1993/miniforge3/envs/alexandros/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/alexsta1993/alexandros/Marmot

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alexsta1993/alexandros/Marmot/build

# Include any dependencies generated for this target.
include CMakeFiles/TestB4.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/TestB4.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/TestB4.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestB4.dir/flags.make

CMakeFiles/TestB4.dir/codegen:
.PHONY : CMakeFiles/TestB4.dir/codegen

CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o: CMakeFiles/TestB4.dir/flags.make
CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o: /home/alexsta1993/alexandros/Marmot/modules/materials/B4/test/test.cpp
CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o: CMakeFiles/TestB4.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/alexsta1993/alexandros/Marmot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o"
	/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o -MF CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o.d -o CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o -c /home/alexsta1993/alexandros/Marmot/modules/materials/B4/test/test.cpp

CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.i"
	/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexsta1993/alexandros/Marmot/modules/materials/B4/test/test.cpp > CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.i

CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.s"
	/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexsta1993/alexandros/Marmot/modules/materials/B4/test/test.cpp -o CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.s

# Object files for target TestB4
TestB4_OBJECTS = \
"CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o"

# External object files for target TestB4
TestB4_EXTERNAL_OBJECTS =

bin/TestB4: CMakeFiles/TestB4.dir/modules/materials/B4/test/test.cpp.o
bin/TestB4: CMakeFiles/TestB4.dir/build.make
bin/TestB4: CMakeFiles/TestB4.dir/compiler_depend.ts
bin/TestB4: lib/libMarmot.so
bin/TestB4: CMakeFiles/TestB4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/alexsta1993/alexandros/Marmot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/TestB4"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestB4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestB4.dir/build: bin/TestB4
.PHONY : CMakeFiles/TestB4.dir/build

CMakeFiles/TestB4.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestB4.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestB4.dir/clean

CMakeFiles/TestB4.dir/depend:
	cd /home/alexsta1993/alexandros/Marmot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alexsta1993/alexandros/Marmot /home/alexsta1993/alexandros/Marmot /home/alexsta1993/alexandros/Marmot/build /home/alexsta1993/alexandros/Marmot/build /home/alexsta1993/alexandros/Marmot/build/CMakeFiles/TestB4.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/TestB4.dir/depend

