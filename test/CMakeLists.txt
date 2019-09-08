find_package(GTest REQUIRED)
add_executable(gtest_converter gtest_converter.cpp)
target_link_libraries(gtest_converter
  PRIVATE
    GTest::GTest
    GTest::Main
)
