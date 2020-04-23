cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 11)
project(SVR)

set(TESTING_BINARY ${CMAKE_PROJECT_NAME}_tests)
set(BENCHMARK_BINARY ${CMAKE_PROJECT_NAME}_benchmarks)
set(TESTING_SOURCE_FILES ../spherical_volume_rendering_util.cpp SVR_tests.cpp)
set(BENCHMARK_SOURCE_FILES ../spherical_volume_rendering_util.cpp SVR_benchmarks.cpp)

add_definitions(-DNDEBUG)
add_subdirectory(googletest)
add_subdirectory(benchmark)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR})

add_executable(${BENCHMARK_BINARY} ${BENCHMARK_SOURCE_FILES})
add_executable(${TESTING_BINARY} ${TESTING_SOURCE_FILES})

target_link_libraries(${BENCHMARK_BINARY} benchmark::benchmark)
target_link_libraries(${TESTING_BINARY} gtest gtest_main gmock gmock_main)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3 -march=native -flto -fno-signed-zeros -funroll-loops -Wall -Wextra")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)