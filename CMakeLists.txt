cmake_minimum_required(VERSION 3.25)
project(Inviscid_Euler_Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Explicitly set Boost_DIR so CMake knows where to find Boost configuration files.
# (Adjust this path to the actual location where vcpkg installs Boost config files.)
set(Boost_DIR "C:/Users/jmmcl/OneDrive/Desktop/Project/C_Plus_Plus_Libraries/vcpkg-master/installed/x64-mingw-static/share/boost")

# Find Boost with the iostreams component using the config mode.
find_package(Boost REQUIRED COMPONENTS iostreams CONFIG)
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost iostreams library: ${Boost_LIBRARIES}")

# Fetch fmt using FetchContent (optional if your project uses fmt).
include(FetchContent)
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 11.1.4
)
FetchContent_MakeAvailable(fmt)

# Create the executable target.
add_executable(Inviscid_Euler_Solver main.cpp)

# Add include directories: the gnuplot-iostream header folder and Boost’s include dirs.
target_include_directories(Inviscid_Euler_Solver PRIVATE
  "C:/Users/jmmcl/OneDrive/Desktop/Project/C_Plus_Plus_Libraries/gnuplot-iostream-master"
  ${Boost_INCLUDE_DIRS}
)

# Link against fmt and Boost Iostreams.
target_link_libraries(Inviscid_Euler_Solver PRIVATE
  fmt::fmt
  Boost::iostreams
)