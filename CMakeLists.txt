cmake_minimum_required(VERSION 3.25)
project(Inviscid_Euler_Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(Matplot++ CONFIG REQUIRED)

# Adds ComputationalGrid as a Library
add_library(ComputationalGridLib
    ComputationalGrid.cpp
)

target_link_libraries(ComputationalGridLib PUBLIC
    fmt::fmt
    Eigen3::Eigen
    Matplot++::matplot
)

add_executable(Inviscid_Euler_Solver main.cpp)

target_link_libraries(Inviscid_Euler_Solver PRIVATE
    ComputationalGridLib
)
