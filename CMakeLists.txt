cmake_minimum_required(VERSION 3.27)
project(DiVerSeS)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#SET(CMAKE_CXX_FLAGS "-std=c++17 -O0") # profiling
#SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_DEBUG -Wreturn-type -fsanitize=address -O0") # debugging
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -march=native")


#include(cmake/DownloadCPSAT.cmake)
include(cmake/cpsat_minimal.cmake)
#include(cmake/deps_gtest_benchmark.cmake)

find_package(OpenMP REQUIRED)

include_directories(
        include
        include/graphs
        include/combinatorics
        include/utils
        include/datastructures
#        include/CONTESTS/PACE22/*
#        include/CONTESTS/PACE22/*/*
        )

file(GLOB SOURCES
        src/graphs/GraphReader.cpp
        src/graphs/GraphWriter.cpp
        src/graphs/GraphUtils.cpp
        src/graphs/GraphInducer.cpp
        src/graphs/graphtraversals/BFS.cpp
        src/graphs/cycles/CycleCounter.cpp

        src/graphs/components/ConnctedComponents.cpp
        src/graphs/cliques/CliqueExtension.cpp
        src/graphs/cliques/CliqueUtils.cpp

        src/graphs/VertexCover/VCUtils.cpp
        src/graphs/VertexCover/kernelization/KernelizerVC.cpp
        src/graphs/matching/MaxMatchBipartite.cpp

        src/graphs/toposort/TopoSort.cpp
        src/graphs/scc/StronglyConnectedComponents.cpp

        src/combinatorics/CombinatoricUtils.cpp

        src/datastructures/RandomSelectionSet.cpp

        src/utils/RandomNumberGenerators.cpp
        src/utils/StandardUtils.cpp
        src/utils/MemoryUtils.cpp
        src/utils/TimeMeasurer.cpp
        src/utils/Stopwatch.cpp

        #PACE22
        src/CONTESTS/PACE22/Config.cpp
        src/CONTESTS/PACE22/heur/AgentFlow.cpp
        src/CONTESTS/PACE22/heur/DFVSSolverH.cpp
        src/CONTESTS/PACE22/heur/VCImprover.cpp
        src/CONTESTS/PACE22/heur/SALS.cpp
        src/CONTESTS/PACE22/heur/SALS2.cpp
        src/CONTESTS/PACE22/heur/SALS3.cpp
        src/CONTESTS/PACE22/heur/SALS4.cpp
        src/CONTESTS/PACE22/heur/SALS5.cpp
        src/CONTESTS/PACE22/heur/DFVSImprover.cpp
        src/CONTESTS/PACE22/heur/HittingSetLS.cpp
        src/CONTESTS/PACE22/Utils.cpp
        src/CONTESTS/PACE22/Reducer.cpp
        src/CONTESTS/PACE22/exact/DFVSSolverE.cpp
        src/CONTESTS/PACE22/exact/CliqueCoverLS.cpp

        src/CONTESTS/PACE22/cpsat-experiment/IntGenerator.h
        src/CONTESTS/PACE22/cpsat-experiment/IntGenerator.cpp
        )


add_executable(DiVerSeS "src/CONTESTS/PACE22/Pace22.cpp" ${SOURCES})
target_link_libraries( DiVerSeS stdc++fs cpsat::cpsat OpenMP::OpenMP_CXX)
#target_link_libraries( DiVerSeS stdc++fs ortools boost_system OpenMP::OpenMP_CXX)



add_executable(CpsatExp1TestGenerator ${SOURCES}
        src/CONTESTS/PACE22/cpsat-experiment/InstancesGenerator.cpp
        src/CONTESTS/PACE22/cpsat-experiment/InstancesGenerator.h
        src/CONTESTS/PACE22/cpsat-experiment/TestGenerator.h
        src/CONTESTS/PACE22/cpsat-experiment/TestGenerator.cpp
)
target_link_libraries( CpsatExp1TestGenerator stdc++fs OpenMP::OpenMP_CXX)

add_executable(CpsatExp1 ${SOURCES}
        src/CONTESTS/PACE22/cpsat-experiment/ExpConfig.cpp
        src/CONTESTS/PACE22/cpsat-experiment/ExpConfig.h
        src/CONTESTS/PACE22/cpsat-experiment/main.cpp
        src/CONTESTS/PACE22/cpsat-experiment/CpsatExp1.cpp
        src/CONTESTS/PACE22/cpsat-experiment/CpsatExp1.h
)
target_link_libraries( CpsatExp1 stdc++fs cpsat::cpsat OpenMP::OpenMP_CXX ortools::ortools)