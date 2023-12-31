cmake_minimum_required(VERSION 3.10)
project(utils)

message(STATUS "[${PROJECT_NAME}] Building....")

add_library(utils SHARED ${PROJECT_SOURCE_DIR}/src/utils.cpp)
target_link_libraries(utils matplotlib_cpp)

add_library(kdtree SHARED ${PROJECT_SOURCE_DIR}/src/KDTree.cpp)

add_executable(ipopt_solve_test ${PROJECT_SOURCE_DIR}/src/ipopt_solve_test.cpp)
target_link_libraries(ipopt_solve_test ipopt)

add_executable(osqp_solve_test ${PROJECT_SOURCE_DIR}/src/osqp_solve_test.cpp)
target_link_libraries(osqp_solve_test OsqpEigen::OsqpEigen)
