cmake_minimum_required(VERSION 3.14)
project(reachSolver)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)

add_subdirectory(src)
#add_subdirectory(tests)
#add_subdirectory(examples)
#add_subdirectory(inner_example)
#add_subdirectory(CORA_example)
add_subdirectory(extra_inner_example)