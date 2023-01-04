cmake_minimum_required(VERSION 3.14)
project(reach_solver_test)

set(CMAKE_CXX_STANDARD 14)

set(TARGET_NAME overElectroOsc)
#test      good
#example2poly
#example2zono
#tank_use_poly
#convert2polygon
#plotter_test
#test_von
#test_memory
#capd_test
#PID_test
#PID_test_simplify
#polyzono_test
#boost_test
#example2poly_compare
#example2poly_border
#example1poly_border
#example2zono_bound
#example1zono_bound
#electroOscZono_bound
#electroOscPoly_bound
#example3
#linear_example1
#glpk_test
#inner_test
#under2016Ex
#under2016Ex8
#underVanderPol
#underTank
#neural_bound_example
#overElectroOsc
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
