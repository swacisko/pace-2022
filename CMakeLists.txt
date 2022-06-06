cmake_minimum_required(VERSION 3.14)
project(DiFVerSeS)

set(CMAKE_CXX_STANDARD 17)


#SET(CMAKE_CXX_FLAGS "-std=c++17 -O0") # profiling
#SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_DEBUG -std=c++17 -Wreturn-type -O0") # debugging

SET(CMAKE_CXX_FLAGS "-std=c++17 -O3")



include_directories(
        include
        include/graphs
        include/combinatorics
        include/algebra
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
        )


add_executable(DiFVerSeS "src/CONTESTS/PACE22/Pace22.cpp" ${SOURCES})
target_link_libraries( DiFVerSeS stdc++fs )