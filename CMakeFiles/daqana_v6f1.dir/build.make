# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/ana/Programs/Daqana_Code/daqana6U

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ana/Programs/Daqana_Code

# Include any dependencies generated for this target.
include CMakeFiles/daqana_v6f1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/daqana_v6f1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/daqana_v6f1.dir/flags.make

CMakeFiles/daqana_v6f1.dir/daqana_v6f1.o: CMakeFiles/daqana_v6f1.dir/flags.make
CMakeFiles/daqana_v6f1.dir/daqana_v6f1.o: daqana6U/daqana_v6f1.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ana/Programs/Daqana_Code/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/daqana_v6f1.dir/daqana_v6f1.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/daqana_v6f1.dir/daqana_v6f1.o -c /home/ana/Programs/Daqana_Code/daqana6U/daqana_v6f1.cc

CMakeFiles/daqana_v6f1.dir/daqana_v6f1.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqana_v6f1.dir/daqana_v6f1.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ana/Programs/Daqana_Code/daqana6U/daqana_v6f1.cc > CMakeFiles/daqana_v6f1.dir/daqana_v6f1.i

CMakeFiles/daqana_v6f1.dir/daqana_v6f1.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqana_v6f1.dir/daqana_v6f1.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ana/Programs/Daqana_Code/daqana6U/daqana_v6f1.cc -o CMakeFiles/daqana_v6f1.dir/daqana_v6f1.s

CMakeFiles/daqana_v6f1.dir/RootFileManager.o: CMakeFiles/daqana_v6f1.dir/flags.make
CMakeFiles/daqana_v6f1.dir/RootFileManager.o: daqana6U/RootFileManager.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ana/Programs/Daqana_Code/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/daqana_v6f1.dir/RootFileManager.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/daqana_v6f1.dir/RootFileManager.o -c /home/ana/Programs/Daqana_Code/daqana6U/RootFileManager.cc

CMakeFiles/daqana_v6f1.dir/RootFileManager.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqana_v6f1.dir/RootFileManager.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ana/Programs/Daqana_Code/daqana6U/RootFileManager.cc > CMakeFiles/daqana_v6f1.dir/RootFileManager.i

CMakeFiles/daqana_v6f1.dir/RootFileManager.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqana_v6f1.dir/RootFileManager.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ana/Programs/Daqana_Code/daqana6U/RootFileManager.cc -o CMakeFiles/daqana_v6f1.dir/RootFileManager.s

# Object files for target daqana_v6f1
daqana_v6f1_OBJECTS = \
"CMakeFiles/daqana_v6f1.dir/daqana_v6f1.o" \
"CMakeFiles/daqana_v6f1.dir/RootFileManager.o"

# External object files for target daqana_v6f1
daqana_v6f1_EXTERNAL_OBJECTS =

daqana_v6f1: CMakeFiles/daqana_v6f1.dir/daqana_v6f1.o
daqana_v6f1: CMakeFiles/daqana_v6f1.dir/RootFileManager.o
daqana_v6f1: CMakeFiles/daqana_v6f1.dir/build.make
daqana_v6f1: /home/ana/root/lib/libCore.so
daqana_v6f1: /home/ana/root/lib/libImt.so
daqana_v6f1: /home/ana/root/lib/libRIO.so
daqana_v6f1: /home/ana/root/lib/libNet.so
daqana_v6f1: /home/ana/root/lib/libHist.so
daqana_v6f1: /home/ana/root/lib/libGraf.so
daqana_v6f1: /home/ana/root/lib/libGraf3d.so
daqana_v6f1: /home/ana/root/lib/libGpad.so
daqana_v6f1: /home/ana/root/lib/libROOTDataFrame.so
daqana_v6f1: /home/ana/root/lib/libTree.so
daqana_v6f1: /home/ana/root/lib/libTreePlayer.so
daqana_v6f1: /home/ana/root/lib/libRint.so
daqana_v6f1: /home/ana/root/lib/libPostscript.so
daqana_v6f1: /home/ana/root/lib/libMatrix.so
daqana_v6f1: /home/ana/root/lib/libPhysics.so
daqana_v6f1: /home/ana/root/lib/libMathCore.so
daqana_v6f1: /home/ana/root/lib/libThread.so
daqana_v6f1: /home/ana/root/lib/libMultiProc.so
daqana_v6f1: CMakeFiles/daqana_v6f1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ana/Programs/Daqana_Code/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable daqana_v6f1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/daqana_v6f1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/daqana_v6f1.dir/build: daqana_v6f1

.PHONY : CMakeFiles/daqana_v6f1.dir/build

CMakeFiles/daqana_v6f1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/daqana_v6f1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/daqana_v6f1.dir/clean

CMakeFiles/daqana_v6f1.dir/depend:
	cd /home/ana/Programs/Daqana_Code && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ana/Programs/Daqana_Code/daqana6U /home/ana/Programs/Daqana_Code/daqana6U /home/ana/Programs/Daqana_Code /home/ana/Programs/Daqana_Code /home/ana/Programs/Daqana_Code/CMakeFiles/daqana_v6f1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/daqana_v6f1.dir/depend

