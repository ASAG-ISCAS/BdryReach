cmake_minimum_required(VERSION 3.01)

project(reach_solver VERSION 1.0.0)


find_package(Eigen3 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED QUIET)

option(SEPARATE_TEMPLATE "separate template include and implementation" ON)

add_library(${PROJECT_NAME})


set_target_properties(${PROJECT_NAME}
        PROPERTIES
        LINKER_LANGUAGE
        CXX)


target_compile_options(${PROJECT_NAME}
        PUBLIC
        -D__USE_FILIB__)
target_compile_definitions(${PROJECT_NAME}
        PUBLIC
        SEPARATE_TEMPLATE
        _USE_MATH_DEFINES)

target_include_directories(${PROJECT_NAME} PUBLIC
        ${EIGEN3_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
        ${CMAKE_CURRENT_LIST_DIR}/../include)

add_subdirectory(dynamics)
add_subdirectory(zonotope)
add_subdirectory(plotter)

install(TARGETS ${PROJECT_NAME}
        #EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include
        )

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../include/dynamics ${CMAKE_CURRENT_LIST_DIR}/../include/zonotope
        DESTINATION include # target directory
)

#[[install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/dynamics ${CMAKE_CURRENT_LIST_DIR}/zonotope
        DESTINATION src # target directory
)]]

