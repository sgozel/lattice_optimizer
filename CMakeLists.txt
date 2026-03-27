cmake_minimum_required(VERSION 3.21)
project(lattice_optimizer VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g3")
set(COMMON_COMPILE_OPTIONS -Wall)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)

configure_file(data.json data.json)

set(LO_SOURCES
    ${CMAKE_SOURCE_DIR}/src/main.cpp
)

add_executable(main ${LO_SOURCES})
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)
