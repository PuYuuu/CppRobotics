cmake_minimum_required(VERSION 3.10)
project(CppRobotics)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "-std=c++11")

# add_definitions(-DWITHOUT_NUMPY=true)

find_package(PythonLibs REQUIRED)
find_package(fmt REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories( "/usr/include/eigen3")

include_directories(${CMAKE_CURRENT_LIST_DIR})

message("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
message("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

add_executable(move_to_pose move_to_pose/move_to_pose.cpp utils/utils.cpp)
target_link_libraries(move_to_pose fmt::fmt ${PYTHON_LIBRARIES})