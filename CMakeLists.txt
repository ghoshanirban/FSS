cmake_minimum_required(VERSION 3.22)
project(FastSparseSpanner)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fexceptions -Wall -Wextra -Wno-unused  -lgmp -lmpfr")
set(CMAKE_BUILD_TYPE "Release")

add_executable(FastSparseSpanner main.cpp
        PointGenerators.h
        Utilities.h
        FSS/GreedySpanner.h
        FSS/WSPDSpanner.h
        FSS/StretchFactorCalculator.h
        FSS/FSS.h
        measurements/diameter.h)

target_link_libraries(FastSparseSpanner libgmp.a)
