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
include CMakeFiles/TestMenetreyWillam.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/TestMenetreyWillam.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/TestMenetreyWillam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestMenetreyWillam.dir/flags.make

CMakeFiles/TestMenetreyWillam.dir/codegen:
.PHONY : CMakeFiles/TestMenetreyWillam.dir/codegen

CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o: CMakeFiles/TestMenetreyWillam.dir/flags.make
CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o: /home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp
CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o: CMakeFiles/TestMenetreyWillam.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/alexsta1993/alexandros/Marmot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o"
	/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o -MF CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o.d -o CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o -c /home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp

CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.i"
	/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp > CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.i

CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.s"
	/home/alexsta1993/miniforge3/envs/alexandros/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexsta1993/alexandros/Marmot/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp -o CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.s

# Object files for target TestMenetreyWillam
TestMenetreyWillam_OBJECTS = \
"CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o"

# External object files for target TestMenetreyWillam
TestMenetreyWillam_EXTERNAL_OBJECTS =

bin/TestMenetreyWillam: CMakeFiles/TestMenetreyWillam.dir/modules/core/MarmotMechanicsCore/test/TestMenetreyWillam.cpp.o
bin/TestMenetreyWillam: CMakeFiles/TestMenetreyWillam.dir/build.make
bin/TestMenetreyWillam: CMakeFiles/TestMenetreyWillam.dir/compiler_depend.ts
bin/TestMenetreyWillam: lib/libMarmot.so
bin/TestMenetreyWillam: CMakeFiles/TestMenetreyWillam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/alexsta1993/alexandros/Marmot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/TestMenetreyWillam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestMenetreyWillam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestMenetreyWillam.dir/build: bin/TestMenetreyWillam
.PHONY : CMakeFiles/TestMenetreyWillam.dir/build

CMakeFiles/TestMenetreyWillam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestMenetreyWillam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestMenetreyWillam.dir/clean

CMakeFiles/TestMenetreyWillam.dir/depend:
	cd /home/alexsta1993/alexandros/Marmot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alexsta1993/alexandros/Marmot /home/alexsta1993/alexandros/Marmot /home/alexsta1993/alexandros/Marmot/build /home/alexsta1993/alexandros/Marmot/build /home/alexsta1993/alexandros/Marmot/build/CMakeFiles/TestMenetreyWillam.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/TestMenetreyWillam.dir/depend

