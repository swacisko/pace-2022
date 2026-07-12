cmake_minimum_required(VERSION 3.24)
project(DiVerSeS)

set(CMAKE_CXX_STANDARD 23)


#SET(CMAKE_CXX_FLAGS "-std=c++23 -O0") # profiling
#SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_DEBUG -std=c++23 -Wreturn-type -O0") # debugging

SET(CMAKE_CXX_FLAGS "-std=c++23 -O3")


include_directories(
        include
        include/graphs
        include/combinatorics
        include/utils
        include/datastructures
        include/CONTESTS/PACE22/*
        )

file(GLOB SOURCES
        src/graphs/GraphReader.cpp
        src/graphs/GraphWriter.cpp
        src/graphs/GraphUtils.cpp
        src/graphs/GraphInducer.cpp
        src/graphs/components/ConnctedComponents.cpp
        src/graphs/cliques/CliqueExtension.cpp
        src/graphs/cliques/CliqueUtils.cpp
        src/graphs/VertexCover/VCUtils.cpp
        src/graphs/VertexCover/kernelization/KernelizerVC.cpp
        src/graphs/matching/MaxMatchBipartite.cpp
        src/combinatorics/CombinatoricUtils.cpp
        src/utils/RandomNumberGenerators.cpp
        src/utils/StandardUtils.cpp
        src/utils/MemoryUtils.cpp
        src/utils/Stopwatch.cpp

        src/CONTESTS/PACE22/Utils.cpp
        src/CONTESTS/PACE22/Reducer.cpp
        src/CONTESTS/PACE22/exact/DFVSSolverE.cpp
        src/CONTESTS/PACE22/exact/CliqueCoverLS.cpp
)


add_executable(VCReducer ${SOURCES}
        src/CONTESTS/PACE22/EDReducer.cpp
        src/CONTESTS/PACE22/EDReducer.h
        src/CONTESTS/PACE22/EDTests.cpp
        src/CONTESTS/PACE22/IntGenerator.h
        src/CONTESTS/PACE22/IntGenerator.cpp
)

target_compile_features(VCReducer PRIVATE cxx_std_23)

find_package(ortools CONFIG REQUIRED)
find_package(OpenMP)
target_link_libraries( VCReducer stdc++fs ortools::ortools OpenMP::OpenMP_CXX)
