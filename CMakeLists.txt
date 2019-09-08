cmake_minimum_required(VERSION 3.12)
project(cs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -Wextra")

add_subdirectory(test)
enable_testing()
add_test(NAME CPP_TEST_CONVERTER
         COMMAND gtest_converter)
