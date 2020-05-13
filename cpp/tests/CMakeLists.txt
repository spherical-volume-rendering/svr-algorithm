cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 11)
project(SVR)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

#if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#    add_link_options(-Wall -fprofile-arcs -ftest-coverage -fprofile-instr-generate)
#    add_compile_options(-fprofile-instr-generate -fcoverage-mapping)
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    add_compile_options(--coverage)
#elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GCC")
#    add_compile_options(-fprofile-arcs -ftest-coverage)
#endif ()

set(TESTING_BINARY test_${CMAKE_PROJECT_NAME})
set(TESTING_SOURCE_FILES ../spherical_volume_rendering_util.cpp test_SVR.cpp)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
add_executable(${TESTING_BINARY} ${TESTING_SOURCE_FILES})
target_link_libraries(${TESTING_BINARY} gtest_main gmock_main)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "cpp/-Wall -Wextra")

enable_testing()



