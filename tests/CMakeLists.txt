find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(test_matrix
    test_matrix.cpp
)

target_link_libraries(test_matrix
    PRIVATE matrix
    PRIVATE GTest::gtest_main
)

gtest_discover_tests(test_matrix)