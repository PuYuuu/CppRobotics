cmake_minimum_required(VERSION 3.10)
project(CppRobotics)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Configuring Debug build")
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Configuring Release build")
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/../bin)

add_library(matplotlib_cpp INTERFACE)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
target_link_libraries(matplotlib_cpp INTERFACE
  Python3::Python
  Python3::Module
)

find_package(Python3 COMPONENTS NumPy)
if(Python3_NumPy_FOUND)
  target_link_libraries(matplotlib_cpp INTERFACE Python3::NumPy)
else()
  target_compile_definitions(matplotlib_cpp INTERFACE WITHOUT_NUMPY)
endif()

find_package(fmt REQUIRED)
find_package(OsqpEigen REQUIRED)
include_directories("/usr/include/eigen3")
include_directories(${CMAKE_CURRENT_LIST_DIR})

include_directories(
  ${PROJECT_SOURCE_DIR}/3rdparty
  ${PROJECT_SOURCE_DIR}/3rdparty/KDTree
  ${PROJECT_SOURCE_DIR}/utils/include
  ${PROJECT_SOURCE_DIR}/Control/include
  ${PROJECT_SOURCE_DIR}/PathPlanning/include
  ${PROJECT_SOURCE_DIR}/PathTracking/include
  ${PROJECT_SOURCE_DIR}/Perception/include
)

add_library(kdtree SHARED ${PROJECT_SOURCE_DIR}/3rdparty/KDTree/KDTree.cpp)

add_subdirectory(utils)
add_subdirectory(Control)
add_subdirectory(PathPlanning)
add_subdirectory(PathTracking)
add_subdirectory(Perception)
