cmake_minimum_required(VERSION 3.10)
project(PathPlanning)

message(STATUS "[${PROJECT_NAME}] Building....")

add_library(rs_path SHARED
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/reeds_shepp_path.cpp)

add_library(cubic_spline SHARED
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/cubic_spline.cpp)

add_library(graph_search SHARED
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/GraphSearchPlanner.cpp)

add_executable(cubic_spline_path
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/cubic_spline_path.cpp)
add_dependencies(cubic_spline_path utils cubic_spline)
target_link_libraries(cubic_spline_path utils fmt::fmt cubic_spline)

add_executable(bezier_path
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/bezier_path.cpp)
add_dependencies(bezier_path utils)
target_link_libraries(bezier_path utils fmt::fmt)

add_executable(bspline_path
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/bspline_path.cpp)
add_dependencies(bspline_path utils)
target_link_libraries(bspline_path utils fmt::fmt)

add_executable(dubins_path
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/dubins_path.cpp)
add_dependencies(dubins_path utils)
target_link_libraries(dubins_path utils fmt::fmt)

add_executable(frenet_optimal_trajectory
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/frenet_optimal_trajectory.cpp)
add_dependencies(frenet_optimal_trajectory utils cubic_spline)
target_link_libraries(frenet_optimal_trajectory utils fmt::fmt cubic_spline)

add_executable(reeds_shepp_planner
    ${PROJECT_SOURCE_DIR}/src/CurvesGenerator/reeds_shepp_planner.cpp)
add_dependencies(reeds_shepp_planner utils rs_path)
target_link_libraries(reeds_shepp_planner utils fmt::fmt rs_path)

add_executable(breadth_first_search
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/breadth_first_search.cpp)
add_dependencies(breadth_first_search utils graph_search)
target_link_libraries(breadth_first_search utils fmt::fmt graph_search)

add_executable(depth_first_search
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/depth_first_search.cpp)
add_dependencies(depth_first_search utils graph_search)
target_link_libraries(depth_first_search utils fmt::fmt graph_search)

add_executable(dijkstra ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/dijkstra.cpp)
add_dependencies(dijkstra utils graph_search)
target_link_libraries(dijkstra utils fmt::fmt graph_search)

add_executable(astar ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/astar.cpp)
add_dependencies(astar utils graph_search)
target_link_libraries(astar utils fmt::fmt graph_search)

add_executable(astar_bidirectional
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/astar_bidirectional.cpp)
add_dependencies(astar_bidirectional utils graph_search)
target_link_libraries(astar_bidirectional utils fmt::fmt graph_search)

add_executable(hybrid_astar
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/HybridAstar/hybrid_astar.cpp
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/HybridAstar/dynamic_programming_heuristic.cpp)
add_dependencies(hybrid_astar utils rs_path kdtree)
target_link_libraries(hybrid_astar utils fmt::fmt rs_path kdtree)

add_executable(hybrid_astar_with_trailer
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/HybridAstar/hybrid_astar_with_trailer.cpp
    ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/HybridAstar/dynamic_programming_heuristic.cpp)
add_dependencies(hybrid_astar_with_trailer utils rs_path kdtree)
target_link_libraries(hybrid_astar_with_trailer utils fmt::fmt rs_path kdtree)

add_executable(rrt ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/rrt.cpp)
add_dependencies(rrt utils)
target_link_libraries(rrt utils fmt::fmt)

add_executable(prm ${PROJECT_SOURCE_DIR}/src/GlobalPlanner/probabilistic_road_map.cpp)
add_dependencies(prm utils)
target_link_libraries(prm utils fmt::fmt)

add_executable(dynamic_window_approach
    ${PROJECT_SOURCE_DIR}/src/LocalPlanner/dynamic_window_approach.cpp)
add_dependencies(dynamic_window_approach utils)
target_link_libraries(dynamic_window_approach utils fmt::fmt)

add_executable(potential_field_planning
    ${PROJECT_SOURCE_DIR}/src/LocalPlanner/potential_field_planning.cpp)
add_dependencies(potential_field_planning utils)
target_link_libraries(potential_field_planning utils fmt::fmt)

add_executable(lattice_planner
    ${PROJECT_SOURCE_DIR}/src/LocalPlanner/lattice_planner.cpp)
add_dependencies(lattice_planner utils cubic_spline)
target_link_libraries(lattice_planner utils fmt::fmt cubic_spline)
