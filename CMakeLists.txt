cmake_minimum_required(VERSION 3.12)
project(cs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -Wextra")

set(GTEST_INCLUDE $ENV{HOME}/opt/googletest-1.8.0/googletest/include)
set(GTEST_LINK $ENV{HOME}/opt/googletest-1.8.0/build/googlemock/gtest)
include_directories(${GTEST_INCLUDE})
link_directories(${GTEST_LINK})
add_subdirectory(test)
enable_testing()
add_test(NAME CPP_TEST_CONVERTER
        COMMAND gtest_converter)