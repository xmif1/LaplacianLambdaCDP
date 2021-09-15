cmake_minimum_required(VERSION 3.17)
project(LaplacianLambdaCDP)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include)
add_executable(LaplacianLambdaCDP main.cpp LaplacianLambdaCDP.h)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse2 -DNDEBUG") # specify compiler optimisations for Eigen library