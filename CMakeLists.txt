cmake_minimum_required(VERSION 3.10.2)
set(TARGET "6dof-rigid-body-solver")
project(${TARGET} LANGUAGES CXX DESCRIPTION "6oF-rigid-body-solver")

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(Eigen3 3.2 REQUIRED)

find_package(precice REQUIRED CONFIG)
add_executable(
	${TARGET} 
	main.cpp
	Structure0815.cpp
	Tests.cpp)

target_link_libraries(${TARGET} PRIVATE precice::precice Eigen3::Eigen)

target_include_directories(${TARGET} PRIVATE $ENV{PRECICE_ROOT}/src)
target_include_directories(${TARGET} PRIVATE $ENV{PRECICE_ROOT}/thirdparty/prettyprint/include)
set_target_properties(${TARGET} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED TRUE
    CXX_EXTENSIONS FALSE
    )
