cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 11)
project(svr_tests)
set(BINARY run_${CMAKE_PROJECT_NAME})
set(SOURCE_FILES ../spherical_volume_rendering_util.h ../spherical_volume_rendering_util.cpp spherical_volume_rendering_test.cpp)

add_subdirectory(googletest)
include_directories(${../cpp})
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR})

add_executable(${BINARY} ${SOURCE_FILES})
add_test(NAME ${BINARY} COMMAND ${BINARY})
target_link_libraries(${BINARY} gtest gtest_main gmock gmock_main)