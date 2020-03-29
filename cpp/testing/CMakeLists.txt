cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 11)
project(testing)
set(SOURCE_FILES ../spherical_volume_rendering_util.h ../spherical_volume_rendering_util.cpp spherical_volume_rendering_test.cpp)

add_subdirectory(googletest)
include_directories(${../cpp})
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(spherical_volume_rendering_test ${SOURCE_FILES})
target_link_libraries(spherical_volume_rendering_test gtest gtest_main)