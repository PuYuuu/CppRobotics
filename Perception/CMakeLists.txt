cmake_minimum_required(VERSION 3.10)
project(Perception)

message(STATUS "[${PROJECT_NAME}] Building....")

add_executable(ekf_location
    ${PROJECT_SOURCE_DIR}/src/extend_kalman_filter_location.cpp)
add_dependencies(ekf_location utils)
target_link_libraries(ekf_location utils fmt::fmt)

add_executable(rectangle_fitting ${PROJECT_SOURCE_DIR}/src/rectangle_fitting.cpp)
add_dependencies(rectangle_fitting utils)
target_link_libraries(rectangle_fitting utils fmt::fmt)
