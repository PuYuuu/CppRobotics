cmake_minimum_required(VERSION 3.10)
project(Control)

message(STATUS "[${PROJECT_NAME}] Building....")

add_executable(move_to_pose ${PROJECT_SOURCE_DIR}/src/move_to_pose.cpp)
add_dependencies(move_to_pose utils)
target_link_libraries(move_to_pose utils fmt::fmt)

add_executable(move_to_pose_robot ${PROJECT_SOURCE_DIR}/src/move_to_pose_robot.cpp)
add_dependencies(move_to_pose_robot utils)
target_link_libraries(move_to_pose_robot utils fmt::fmt)
