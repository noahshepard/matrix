#include "matrix/matrix.hpp"
#include "vector/vector.hpp"
#include <gtest/gtest.h>

TEST(VectorOperations, VectorAddition) {
  linalg::vector a({2, 0, 1});
  linalg::vector b({3, 1, 0});

  linalg::vector c({5, 1, 1});

  EXPECT_EQ(a + b, c);
}

TEST(VectorOperations, VectorSubtraction) {
  linalg::vector a((std::vector<double>){2, 0, 1});
  linalg::vector b({3, 1, 0});

  linalg::vector c({-1, -1, 1});

  EXPECT_EQ(a - b, c);
}

TEST(VectorOperations, ScalarMultiplication) {
  linalg::vector a({2, 0, 1});

  linalg::vector b({4, 0, 2});

  linalg::vector c({1, 0, 0.5});

  linalg::vector d({-2, 0, -1});

  EXPECT_TRUE(a * 2 == b);
  EXPECT_TRUE(0.5 * a == c);
  EXPECT_TRUE(a / -1 == d);
}

TEST(VectorMultiplication, RotationMat) {
  linalg::matrix rot({{1, 0, 0}, {0, 0, -1}, {0, 1, 0}});

  linalg::vector a({0, 1, 0});
  linalg::vector a1({0, 0, 1});

  linalg::vector b({1, 1, 1});
  linalg::vector b1({1, -1, 1});

  EXPECT_EQ(rot * a, a1);
  EXPECT_EQ(rot * b, b1);
}