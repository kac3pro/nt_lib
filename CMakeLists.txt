cmake_minimum_required(VERSION 3.16)
project(nt_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_subdirectory(googletest)

add_library(nt INTERFACE)

# Make your current source‐dir visible to anything that links nt:
target_include_directories(nt
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Tell CMake that nt.hpp is a “source” of this INTERFACE target
target_sources(nt
  INTERFACE
    nt.hpp
)

add_executable(test_runner test.cpp)

target_link_libraries(test_runner
  PRIVATE
    nt                 # header-only library
    GTest::gtest_main  # from add_subdirectory(googletest)
)

include(GoogleTest)
gtest_discover_tests(test_runner)