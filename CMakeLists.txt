cmake_minimum_required(VERSION 3.10.2)
set(TARGET "6dof-rigid-body-solver")
project(${TARGET} LANGUAGES CXX DESCRIPTION "6oF-rigid-body-solver")


find_package(precice REQUIRED CONFIG)
add_executable(
	${TARGET} 
	main.cpp
	Structure0815.cpp
	Tests.cpp)

target_link_libraries(${TARGET} PRIVATE precice::precice)

target_include_directories(${TARGET} PRIVATE $ENV{PRECICE_ROOT}/src)
target_include_directories(${TARGET} PRIVATE $ENV{PRECICE_ROOT}/thirdparty/prettyprint/include)
set_target_properties(${TARGET} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED TRUE
    CXX_EXTENSIONS FALSE
    )
