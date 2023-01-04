cmake_minimum_required(VERSION 3.14)
project(reach_solver_test)

set(CMAKE_CXX_STANDARD 14)

set(TARGET_NAME ElectroOsc)
#Brusselator
#ElectroOsc

find_package(Eigen3 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED QUIET)

add_executable(${TARGET_NAME} ${TARGET_NAME}.cpp)
target_link_libraries(${TARGET_NAME} PUBLIC
        #prim
        ${PYTHON_LIBRARIES}
        /usr/local/lib/libcapd.so
        #        capd
        Eigen3::Eigen
        glpk
        reach_solver)
