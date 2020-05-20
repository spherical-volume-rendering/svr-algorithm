cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 11)
project(SVR)

add_definitions(-DNDEBUG) # Run benchmarks in release mode.

include(FetchContent)
FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG master)
FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG master)
FetchContent_MakeAvailable(
        googletest
        googlebenchmark)

set(BENCHMARK_BINARY benchmark_${CMAKE_PROJECT_NAME})
set(BENCHMARK_SOURCE_FILES ../spherical_volume_rendering_util.cpp benchmark_svr.cpp)

add_executable(${BENCHMARK_BINARY} ${BENCHMARK_SOURCE_FILES})

target_link_libraries(${BENCHMARK_BINARY} benchmark::benchmark)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3 -march=native -flto -fno-signed-zeros -funroll-loops -Wall -Wextra")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)