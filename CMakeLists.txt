cmake_minimum_required(VERSION 3.27)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(fastsim)

list(APPEND CMAKE_PREFIX_PATH "/home/chlorinepentoxide/Qt/6.6.1/gcc_64")

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        cxxplot
        GIT_REPOSITORY https://github.com/USNavalResearchLaboratory/cxxplot.git
        GIT_TAG        v0.4.1
)

FetchContent_MakeAvailable( cxxplot )

add_executable(fastsim fastsim2024.cpp experiment_1.cpp experiment_2.cpp
                        lib/jacobian.h lib/range.h lib/types.h
                        solvers/euler.h solvers/newton_raphson.h solvers/runge_kutta.h)
target_link_libraries(fastsim PRIVATE cxxplot::cxxplot)

