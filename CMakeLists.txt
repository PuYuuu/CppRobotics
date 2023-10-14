cmake_minimum_required(VERSION 3.10)
project(CppRobotics)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++17")

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

add_executable(breadth_first_search PathPlanning/GlobalPlanner/breadth_first_search.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(breadth_first_search fmt::fmt ${PYTHON_LIBRARIES})

add_executable(depth_first_search PathPlanning/GlobalPlanner/depth_first_search.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(depth_first_search fmt::fmt ${PYTHON_LIBRARIES})

add_executable(dijkstra PathPlanning/GlobalPlanner/dijkstra.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(dijkstra fmt::fmt ${PYTHON_LIBRARIES})

add_executable(astar PathPlanning/GlobalPlanner/astar.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(astar fmt::fmt ${PYTHON_LIBRARIES})

add_executable(astar_bidirectional PathPlanning/GlobalPlanner/astar_bidirectional.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(astar_bidirectional fmt::fmt ${PYTHON_LIBRARIES})

add_executable(rrt PathPlanning/GlobalPlanner/rrt.cpp utils/utils.cpp)
target_link_libraries(rrt fmt::fmt ${PYTHON_LIBRARIES})

add_executable(dynamic_window_approach 
    PathPlanning/LocalPlanner/dynamic_window_approach.cpp utils/utils.cpp)
target_link_libraries(dynamic_window_approach fmt::fmt ${PYTHON_LIBRARIES})

add_executable(pure_pursuit PathTracking/pure_pursuit.cpp utils/utils.cpp)
target_link_libraries(pure_pursuit fmt::fmt ${PYTHON_LIBRARIES})
