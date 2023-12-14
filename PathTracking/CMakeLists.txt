cmake_minimum_required(VERSION 3.10)
project(PathTracking)

message(STATUS "[${PROJECT_NAME}] Building....")

add_executable(pure_pursuit
    ${PROJECT_SOURCE_DIR}/src/pure_pursuit.cpp)
add_dependencies(pure_pursuit utils)
target_link_libraries(pure_pursuit utils fmt::fmt)

add_executable(stanley_controller
    ${PROJECT_SOURCE_DIR}/src/stanley_controller.cpp)
add_dependencies(stanley_controller utils cubic_spline)
target_link_libraries(stanley_controller utils fmt::fmt cubic_spline)

add_executable(lqr_with_cartesian
    ${PROJECT_SOURCE_DIR}/src/lqr_with_cartesian.cpp)
add_dependencies(lqr_with_cartesian utils cubic_spline)
target_link_libraries(lqr_with_cartesian utils fmt::fmt cubic_spline)

add_executable(lqr_with_frenet
    ${PROJECT_SOURCE_DIR}/src/lqr_with_frenet.cpp)
add_dependencies(lqr_with_frenet utils cubic_spline)
target_link_libraries(lqr_with_frenet utils fmt::fmt cubic_spline)

add_executable(model_predictive_control
    ${PROJECT_SOURCE_DIR}/src/model_predictive_control.cpp)
add_dependencies(model_predictive_control utils cubic_spline)
target_link_libraries(model_predictive_control utils fmt::fmt cubic_spline ipopt)
