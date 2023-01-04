cmake_minimum_required(VERSION 3.14)
project(reach_solver_test)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED QUIET)

set(LINKING_LIBRARIES
        #prim
        ${PYTHON_LIBRARIES}
        /usr/local/lib/libcapd.so
        #        capd
        Eigen3::Eigen
        glpk
        reach_solver)

add_executable(over2Rrobot over2Rrobot.cpp)
target_link_libraries(over2Rrobot PUBLIC
        ${LINKING_LIBRARIES})

# add_executable(overTank overTank.cpp)
# target_link_libraries(overTank PUBLIC
#         ${LINKING_LIBRARIES})

# add_executable(overVanderPol overVanderPol.cpp)
# target_link_libraries(overVanderPol 
#         PUBLIC
#         ${LINKING_LIBRARIES})

# add_executable(overElectroOsc overElectroOsc.cpp)
# target_link_libraries(overElectroOsc 
#         PUBLIC
#         ${LINKING_LIBRARIES})

# add_executable(underVanderPol underVanderPol.cpp)
# target_link_libraries(underVanderPol 
#         PUBLIC
#         ${LINKING_LIBRARIES})

# add_executable(underTank underTank.cpp)
# target_link_libraries(underTank 
#         PUBLIC
#         ${LINKING_LIBRARIES})


# add_executable(underElectroOsc underElectroOsc.cpp)
# target_link_libraries(underElectroOsc 
#         PUBLIC
#         ${LINKING_LIBRARIES})

# add_executable(underBiological underBiological.cpp)
# target_link_libraries(underBiological 
#         PUBLIC
#         ${LINKING_LIBRARIES})
