cmake_minimum_required(VERSION 3.10)
project(CppRobotics)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

find_package(PythonLibs REQUIRED)
find_package(fmt REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories( "/usr/include/eigen3")

include_directories(${CMAKE_CURRENT_LIST_DIR})

add_executable(move_to_pose Control/move_to_pose/move_to_pose.cpp utils/utils.cpp)
target_link_libraries(move_to_pose fmt::fmt ${PYTHON_LIBRARIES})

add_executable(move_to_pose_robot Control/move_to_pose/move_to_pose_robot.cpp utils/utils.cpp)
target_link_libraries(move_to_pose_robot fmt::fmt ${PYTHON_LIBRARIES})

add_executable(breadth_first_search PathPlanning/GraphBased/breadth_first_search.cpp
    PathPlanning/GraphBased/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(breadth_first_search fmt::fmt ${PYTHON_LIBRARIES})

add_executable(depth_first_search PathPlanning/GraphBased/depth_first_search.cpp
    PathPlanning/GraphBased/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(depth_first_search fmt::fmt ${PYTHON_LIBRARIES})

add_executable(dijkstra PathPlanning/GraphBased/dijkstra.cpp
    PathPlanning/GraphBased/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(dijkstra fmt::fmt ${PYTHON_LIBRARIES})

add_executable(astar PathPlanning/GraphBased/astar.cpp
    PathPlanning/GraphBased/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(astar fmt::fmt ${PYTHON_LIBRARIES})
