# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 4.0

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\build

# Include any dependencies generated for this target.
include CMakeFiles/InitializeLib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/InitializeLib.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/InitializeLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/InitializeLib.dir/flags.make

CMakeFiles/InitializeLib.dir/codegen:
.PHONY : CMakeFiles/InitializeLib.dir/codegen

CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj: CMakeFiles/InitializeLib.dir/flags.make
CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj: CMakeFiles/InitializeLib.dir/includes_CXX.rsp
CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj: C:/Users/jmmcl/OneDrive/Desktop/Project/CFD/Inviscid_Euler_Solver/src/Initialize.cpp
CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj: CMakeFiles/InitializeLib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj -MF CMakeFiles\InitializeLib.dir\src\Initialize.cpp.obj.d -o CMakeFiles\InitializeLib.dir\src\Initialize.cpp.obj -c C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\src\Initialize.cpp

CMakeFiles/InitializeLib.dir/src/Initialize.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/InitializeLib.dir/src/Initialize.cpp.i"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\src\Initialize.cpp > CMakeFiles\InitializeLib.dir\src\Initialize.cpp.i

CMakeFiles/InitializeLib.dir/src/Initialize.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/InitializeLib.dir/src/Initialize.cpp.s"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\src\Initialize.cpp -o CMakeFiles\InitializeLib.dir\src\Initialize.cpp.s

# Object files for target InitializeLib
InitializeLib_OBJECTS = \
"CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj"

# External object files for target InitializeLib
InitializeLib_EXTERNAL_OBJECTS =

libInitializeLib.a: CMakeFiles/InitializeLib.dir/src/Initialize.cpp.obj
libInitializeLib.a: CMakeFiles/InitializeLib.dir/build.make
libInitializeLib.a: CMakeFiles/InitializeLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libInitializeLib.a"
	$(CMAKE_COMMAND) -P CMakeFiles\InitializeLib.dir\cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\InitializeLib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/InitializeLib.dir/build: libInitializeLib.a
.PHONY : CMakeFiles/InitializeLib.dir/build

CMakeFiles/InitializeLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\InitializeLib.dir\cmake_clean.cmake
.PHONY : CMakeFiles/InitializeLib.dir/clean

CMakeFiles/InitializeLib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\build C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\build C:\Users\jmmcl\OneDrive\Desktop\Project\CFD\Inviscid_Euler_Solver\build\CMakeFiles\InitializeLib.dir\DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/InitializeLib.dir/depend

