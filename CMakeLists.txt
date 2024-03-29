cmake_minimum_required(VERSION 3.27)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(fastsim)

list(APPEND CMAKE_PREFIX_PATH "/opt/Qt/6.6.1/gcc_64/")

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        cxxplot
        GIT_REPOSITORY https://github.com/USNavalResearchLaboratory/cxxplot.git
        GIT_TAG 0ffe21b123ff86c6a3aac80e2bf5ab6412169b78
)

FetchContent_MakeAvailable( cxxplot )

add_executable(fastsim fastsim2024.cpp experiments/internal/experiment_1.cpp experiments/internal/experiment_2.cpp
                        lib/jacobian.h lib/range.h lib/types.h
                        solvers/euler.h solvers/newton_raphson.h solvers/runge_kutta.h
        experiments/internal/experiment_3.cpp)

add_executable(vanderpol examples/vanderpol.cpp)
add_executable(lorenz examples/lorenz.cpp)
add_executable(utsav_model examples/utsav_model.cpp
        rocket/models/urm_model.h
)
add_executable(liu_model rocket/models/liu_model.cpp)
add_executable(exp_rlcp231 experiments/researchlab/exp_rlcp231.cpp
        rocket/models/lynx/lynx.h)
add_executable(exp_rlcp232 experiments/researchlab/exp_rlcp232.cpp)
add_executable(exp_rlcp233 experiments/researchlab/exp_rlcp233.cpp)
add_executable(exp_rlcp234 experiments/researchlab/exp_rlcp234.cpp)
add_executable(lynx
        rocket/models/lynx/lynx.h
        rocket/models/lynx/lynx_types.h
        rocket/models/lynx/lynx.cpp)

target_link_libraries(fastsim PRIVATE cxxplot::cxxplot)
target_link_libraries(vanderpol PRIVATE cxxplot::cxxplot)
target_link_libraries(lorenz PRIVATE cxxplot::cxxplot)
target_link_libraries(utsav_model PRIVATE cxxplot::cxxplot)
target_link_libraries(liu_model PRIVATE cxxplot::cxxplot)
target_link_libraries(exp_rlcp231 PRIVATE cxxplot::cxxplot)
target_link_libraries(exp_rlcp232 PRIVATE cxxplot::cxxplot)
target_link_libraries(exp_rlcp233 PRIVATE cxxplot::cxxplot)
target_link_libraries(exp_rlcp234 PRIVATE cxxplot::cxxplot)
target_link_libraries(lynx PRIVATE cxxplot::cxxplot)
