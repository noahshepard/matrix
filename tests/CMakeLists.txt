find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(test_linalg
    test_matrix.cpp
    test_vector.cpp
)

target_link_libraries(test_linalg
    PRIVATE linalg
    PRIVATE GTest::gtest_main
)

gtest_discover_tests(test_linalg)