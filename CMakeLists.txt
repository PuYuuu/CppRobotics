cmake_minimum_required(VERSION 3.10)
project(CppRobotics)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "-std=c++17")

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
include_directories( "/usr/include/eigen3")
include_directories(${CMAKE_CURRENT_LIST_DIR})

add_library(RS_Path_LIB SHARED
    ${PROJECT_SOURCE_DIR}/PathPlanning/CurvesGenerator/reeds_shepp_path.cpp)

add_executable(move_to_pose Control/move_to_pose.cpp utils/utils.cpp)
target_link_libraries(move_to_pose fmt::fmt matplotlib_cpp)

add_executable(move_to_pose_robot Control/move_to_pose_robot.cpp utils/utils.cpp)
target_link_libraries(move_to_pose_robot fmt::fmt matplotlib_cpp)

add_executable(cubic_spline_path PathPlanning/CurvesGenerator/cubic_spline_path.cpp 
    PathPlanning/CurvesGenerator/cubic_spline.cpp utils/utils.cpp)
target_link_libraries(cubic_spline_path fmt::fmt matplotlib_cpp)

add_executable(bezier_path PathPlanning/CurvesGenerator/bezier_path.cpp utils/utils.cpp)
target_link_libraries(bezier_path fmt::fmt matplotlib_cpp)

add_executable(bspline_path PathPlanning/CurvesGenerator/bspline_path.cpp utils/utils.cpp)
target_link_libraries(bspline_path fmt::fmt matplotlib_cpp)

add_executable(dubins_path PathPlanning/CurvesGenerator/dubins_path.cpp utils/utils.cpp)
target_link_libraries(dubins_path fmt::fmt matplotlib_cpp)

add_executable(reeds_shepp_planner PathPlanning/CurvesGenerator/reeds_shepp_planner.cpp utils/utils.cpp)
target_link_libraries(reeds_shepp_planner fmt::fmt matplotlib_cpp RS_Path_LIB)

add_executable(breadth_first_search PathPlanning/GlobalPlanner/breadth_first_search.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(breadth_first_search fmt::fmt matplotlib_cpp)

add_executable(depth_first_search PathPlanning/GlobalPlanner/depth_first_search.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(depth_first_search fmt::fmt matplotlib_cpp)

add_executable(dijkstra PathPlanning/GlobalPlanner/dijkstra.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(dijkstra fmt::fmt matplotlib_cpp)

add_executable(hybrid_astar PathPlanning/GlobalPlanner/HybridAstar/hybrid_astar.cpp
    PathPlanning/GlobalPlanner/HybridAstar/dynamic_programming_heuristic.cpp
    utils/KDTree.cpp utils/utils.cpp)
target_link_libraries(hybrid_astar fmt::fmt matplotlib_cpp RS_Path_LIB)

add_executable(hybrid_astar_with_trailer 
    PathPlanning/GlobalPlanner/HybridAstar/hybrid_astar_with_trailer.cpp
    PathPlanning/GlobalPlanner/HybridAstar/dynamic_programming_heuristic.cpp
    utils/KDTree.cpp utils/utils.cpp)
target_link_libraries(hybrid_astar_with_trailer fmt::fmt matplotlib_cpp RS_Path_LIB)

add_executable(astar PathPlanning/GlobalPlanner/astar.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(astar fmt::fmt matplotlib_cpp)

add_executable(astar_bidirectional PathPlanning/GlobalPlanner/astar_bidirectional.cpp
    PathPlanning/GlobalPlanner/GraphSearchPlanner.cpp utils/utils.cpp)
target_link_libraries(astar_bidirectional fmt::fmt matplotlib_cpp)

add_executable(rrt PathPlanning/GlobalPlanner/rrt.cpp utils/utils.cpp)
target_link_libraries(rrt fmt::fmt matplotlib_cpp)

add_executable(prm PathPlanning/GlobalPlanner/probabilistic_road_map.cpp utils/utils.cpp)
target_link_libraries(prm fmt::fmt matplotlib_cpp)

add_executable(dynamic_window_approach 
    PathPlanning/LocalPlanner/dynamic_window_approach.cpp utils/utils.cpp)
target_link_libraries(dynamic_window_approach fmt::fmt matplotlib_cpp)

add_executable(potential_field_planning 
    PathPlanning/LocalPlanner/potential_field_planning.cpp utils/utils.cpp)
target_link_libraries(potential_field_planning fmt::fmt matplotlib_cpp)

add_executable(pure_pursuit PathTracking/pure_pursuit.cpp utils/utils.cpp)
target_link_libraries(pure_pursuit fmt::fmt matplotlib_cpp)

add_executable(stanley_controller PathTracking/stanley_controller.cpp 
    PathPlanning/CurvesGenerator/cubic_spline.cpp utils/utils.cpp)
target_link_libraries(stanley_controller fmt::fmt matplotlib_cpp)

add_executable(ekf_location Perception/extend_kalman_filter_location.cpp utils/utils.cpp)
target_link_libraries(ekf_location fmt::fmt matplotlib_cpp)

add_executable(rectangle_fitting Perception/rectangle_fitting utils/utils.cpp)
target_link_libraries(rectangle_fitting fmt::fmt matplotlib_cpp)
