cmake_minimum_required(VERSION 3.25)
project(Inviscid_Euler_Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# If Boost iostreams was built with compression support, you need zlib & bzip2.
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

# Find Boost with iostreams
find_package(Boost REQUIRED COMPONENTS iostreams)

# Find fmt
find_package(fmt REQUIRED)

# Find Eigen3 - Eigen is header-only, so no linking is necessary.
# Using NO_MODULE tells CMake to look for the Eigen3Config.cmake file.
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Create the executable
add_executable(Inviscid_Euler_Solver main.cpp)

# If using gnuplot-iostream (header-only), just include its directory:
target_include_directories(Inviscid_Euler_Solver PRIVATE
    "C:/Users/jmmcl/OneDrive/Desktop/Project/C_Plus_Plus_Libraries/gnuplot-iostream-master"
)

# Link the libraries found by vcpkg (Boost, ZLIB, BZip2, fmt) and include Eigen3.
target_link_libraries(Inviscid_Euler_Solver PRIVATE
    Boost::iostreams
    ZLIB::ZLIB
    BZip2::BZip2
    fmt::fmt
    Eigen3::Eigen
)
