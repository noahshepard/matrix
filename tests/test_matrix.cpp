#include <gtest/gtest.h>
#include "matrix/matrix.hpp"

constexpr double EPS = 1e-9;

// Helper to compare doubles with EPS tolerance
void EXPECT_DOUBLE_EQ_EPS(double a, double b)
{
    EXPECT_NEAR(a, b, EPS);
}

// Test 1: Simple 2x3 system (unique solution)
TEST(MatrixRREF, SimpleAugmented2x3)
{
    Matrix A({{1, 2, 1, 9},
              {2, -1, 1, 8}});

    A.rref();

    // Expected RREF:
    // [1 0 1/5 5]
    // [0 1 2/5 2]
    EXPECT_DOUBLE_EQ_EPS(A(0, 0), 1.0);
    EXPECT_DOUBLE_EQ_EPS(A(0, 1), 0.0);
    EXPECT_DOUBLE_EQ_EPS(A(1, 0), 0.0);
    EXPECT_DOUBLE_EQ_EPS(A(1, 1), 1.0);
}

// Test 2: 3x4 system with zero row
TEST(MatrixRREF, AugmentedWithZeroRow)
{
    Matrix A({{1, 2, 1, 9},
              {2, 4, 2, 18},
              {0, 0, 0, 0}});

    A.rref();

    // The first two rows should be normalized, last row all zeros
    EXPECT_DOUBLE_EQ_EPS(A(0, 0), 1.0);
    for (size_t r = 1; r < A.rows(); ++r)
        for (size_t c = 0; c < A.cols(); ++c)
            EXPECT_DOUBLE_EQ_EPS(A(r, c), 0.0);
}

// Test 3: Rank-deficient system
TEST(MatrixRREF, RankDeficientAugmented)
{
    Matrix A({{1, 2, 3, 6},
              {2, 4, 6, 12},
              {3, 6, 9, 18}});

    A.rref();

    // Only first row should have pivot = 1, other rows all zeros
    EXPECT_DOUBLE_EQ_EPS(A(0, 0), 1.0);
    for (size_t r = 1; r < A.rows(); ++r)
        for (size_t c = 0; c < A.cols(); ++c)
            EXPECT_DOUBLE_EQ_EPS(A(r, c), 0.0);
}

// Test 4: 3x4 system with unique solution
TEST(MatrixRREF, FullAugmented3x4)
{
    Matrix A({{2, 1, -1, 8},
              {-3, -1, 2, -11},
              {-2, 1, 2, -3}});

    A.rref();

    // Expected solution: x=2, y=3, z=-1
    EXPECT_DOUBLE_EQ_EPS(A(0, 3), 2.0);  // x
    EXPECT_DOUBLE_EQ_EPS(A(1, 3), 3.0);  // y
    EXPECT_DOUBLE_EQ_EPS(A(2, 3), -1.0); // z
}