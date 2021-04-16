cmake_minimum_required (VERSION 3.0)

project(ZGraph)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Numa REQUIRED)

# version number
set (GraphZ_VERSION_MAJOR 0)
set (GraphZ_VERSION_MINOR 1)

# setup compiler
set (CMAKE_C_COMPILER "mpicc")
set (CMAKE_CXX_COMPILER "mpicxx")

# cxx flags
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "-O3")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D PRINT_DEBUG_MSG") # delete this line to prevent printing debug msg

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories(NUMA_INCLUDE_DIR)

add_library(Bitmap src/Bitmap.cpp)
add_library(Random src/Random.cpp)
add_library(Time src/Time.cpp)
add_library(Debug src/Debug.cpp)
add_library(NetworkInterface src/NetworkInterface.cpp)
add_library(Utilities src/Utilities.cpp)
add_library(Graph src/Graph.cpp)
add_library(SharedMemSys src/SharedMemSys.cpp)
add_library(ZGraphInstance src/ZGraphInstance.cpp)
add_library(ErrorLatencyProfile src/ErrorLatencyProfile.cpp)
add_library(Regression src/Regression.cpp)
add_library(GraphFilters src/GraphFilters.cpp)
add_library(SamplerGenerator src/SamplerGenerator.cpp)

# setup dependencies
set (dependencies SamplerGenerator)
set (dependencies ${dependencies} ErrorLatencyProfile)
set (dependencies ${dependencies} ZGraphInstance)
set (dependencies ${dependencies} GraphFilters)
set (dependencies ${dependencies} Graph)
set (dependencies ${dependencies} Bitmap)
set (dependencies ${dependencies} Time)
set (dependencies ${dependencies} SharedMemSys)
set (dependencies ${dependencies} Debug)
set (dependencies ${dependencies} Random)
set (dependencies ${dependencies} NetworkInterface)
set (dependencies ${dependencies} Utilities)
set (dependencies ${dependencies} Regression)
set (dependencies ${dependencies} ${NUMA_LIBRARY})

# tests

add_executable(TestRandom tests/TestRandom.cpp)
target_link_libraries(TestRandom ${dependencies})

add_executable(TestGraph tests/TestGraph.cpp)
target_link_libraries(TestGraph ${dependencies})

add_executable(TestSubgraph tests/TestSubgraph.cpp)
target_link_libraries(TestSubgraph ${dependencies})

add_executable(TestRegression tests/TestRegression.cpp)
target_link_libraries(TestRegression ${dependencies})

add_executable(TestGraphFilter tests/TestGraphFilter.cpp)
target_link_libraries(TestGraphFilter ${dependencies})

add_executable(TestPattern tests/TestPattern.cpp)
target_link_libraries(TestPattern ${dependencies})

add_executable(TestPattern2 tests/TestPattern2.cpp)
target_link_libraries(TestPattern2 ${dependencies})

# applications

add_executable(ThreeChain applications/ThreeChain.cpp)
target_link_libraries(ThreeChain ${dependencies})

add_executable(Triangle applications/Triangle.cpp)
target_link_libraries(Triangle ${dependencies})

#add_executable(FourClique applications/FourClique.cpp)
#target_link_libraries(FourClique ${dependencies})

add_executable(AllPredicateWeightedTriangle applications/AllPredicateWeightedTriangle.cpp)
target_link_libraries(AllPredicateWeightedTriangle ${dependencies})

add_executable(AtLeastOnePredicateWeightedTriangle applications/AtLeastOnePredicateWeightedTriangle.cpp)
target_link_libraries(AtLeastOnePredicateWeightedTriangle ${dependencies})

add_executable(ThreeMotif applications/ThreeMotif.cpp)
target_link_libraries(ThreeMotif ${dependencies})

add_executable(ThreeMotifWithoutCache applications/ThreeMotifWithoutCache.cpp)
target_link_libraries(ThreeMotifWithoutCache ${dependencies})

add_executable(FourChain applications/FourChain.cpp)
target_link_libraries(FourChain ${dependencies})

add_executable(ChainMining applications/ChainMining.cpp)
target_link_libraries(ChainMining ${dependencies})

add_executable(Triangle2 applications/Triangle2.cpp)
target_link_libraries(Triangle2 ${dependencies})

add_executable(FiveStar applications/FiveStar.cpp)
target_link_libraries(FiveStar ${dependencies})

add_executable(FiveHouse applications/FiveHouse.cpp)
target_link_libraries(FiveHouse ${dependencies})
