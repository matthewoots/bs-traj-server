cmake_minimum_required(VERSION 3.0.2)
project(bs_traj_server)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")

file(GLOB INCLUDE_FILES "include/*.h")
file(GLOB SRC_FILES "src/*.cpp")

find_package(Eigen3 REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIRS} 
    libbspline/cpp/include
)

add_subdirectory(libbspline/cpp)

add_executable(${PROJECT_NAME}_test_random_points_node 
    src/test_random_points.cpp
    src/trajectory_server.cpp
)

add_executable(${PROJECT_NAME}_test_concatenate_node 
    src/test_concatenate.cpp
    src/trajectory_server.cpp
)

add_library(${PROJECT_NAME} STATIC ${SRC_FILES} ${INCLUDE_FILES})
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)