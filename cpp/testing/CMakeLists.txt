cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 11)
project(SVR)

set(TESTING_BINARY ${CMAKE_PROJECT_NAME}_tests)
set(BENCHMARK_BINARY ${CMAKE_PROJECT_NAME}_benchmarks)
set(TESTING_SOURCE_FILES ../spherical_volume_rendering_util.cpp SVR_tests.cpp)
set(BENCHMARK_SOURCE_FILES ../spherical_volume_rendering_util.cpp SVR_benchmarks.cpp)

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)

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
