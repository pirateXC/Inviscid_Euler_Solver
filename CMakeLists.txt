cmake_minimum_required(VERSION 3.25)
project(Inviscid_Euler_Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
 # set(CMAKE_BUILD_TYPE Debug CACHE STRING
 #     "Choose the type of build (Debug or Release)" FORCE)
#endif()

# Build‐type‐specific flags
#    - Debug: debug symbols, no optimization
#    - Release: optimize
#set(DEBUG_FLAGS   "-g;-O0")
#set(RELEASE_FLAGS "-O3")

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(Matplot++ CONFIG REQUIRED)

add_library(cfdsolver
    src/GridHandler.cpp
    src/Initialize.cpp
    src/StegerWarming.cpp
    src/PostProcess.cpp)

target_include_directories(cfdsolver
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(cfdsolver PRIVATE
  $<$<CONFIG:Debug>:${DEBUG_FLAGS}>
  $<$<CONFIG:Release>:${RELEASE_FLAGS}>
)

target_link_libraries(cfdsolver
  PUBLIC
    Eigen3::Eigen
    fmt::fmt
    Matplot++::matplot
)

add_executable(Inviscid_Euler_Solver
    src/main.cpp
)

target_compile_options(Inviscid_Euler_Solver PRIVATE
  $<$<CONFIG:Debug>:${DEBUG_FLAGS}>
  $<$<CONFIG:Release>:${RELEASE_FLAGS}>
)

target_link_libraries(Inviscid_Euler_Solver
  PRIVATE
    cfdsolver
)
