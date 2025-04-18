cmake_minimum_required(VERSION 3.25)
project(Inviscid_Euler_Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(Matplot++ CONFIG REQUIRED)

# Include directories for header files
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add GridHandlerLib as a Library 
add_library(GridHandlerLib
    src/GridHandler.cpp
)

# Add InitializeLib as a Library 
add_library(InitializeLib
    src/Initialize.cpp
)

target_link_libraries(GridHandlerLib PUBLIC
    fmt::fmt
    Eigen3::Eigen
    Matplot++::matplot
)

target_link_libraries(InitializeLib PUBLIC
    fmt::fmt
    Eigen3::Eigen
    Matplot++::matplot
)


# Create the main executable, now sourcing main.cpp from the src/ folder
add_executable(Inviscid_Euler_Solver src/main.cpp)

target_link_libraries(Inviscid_Euler_Solver PRIVATE
    GridHandlerLib
    InitializeLib
)
