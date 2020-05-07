cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 11)
project(SVR)

set(TESTING_BINARY test_${CMAKE_PROJECT_NAME})
set(TESTING_SOURCE_FILES ../spherical_volume_rendering_util.cpp test_SVR.cpp)

add_subdirectory(googletest)

add_executable(${TESTING_BINARY} ${TESTING_SOURCE_FILES})

target_link_libraries(${TESTING_BINARY} gtest gtest_main gmock gmock_main)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

