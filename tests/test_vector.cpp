#include <gtest/gtest.h>

#include "matrix/matrix.hpp"
#include "vector/vector.hpp"

TEST(VectorOperations, VectorAddition) {
    linalg::vector a({2, 0, 1});
    linalg::vector b({3, 1, 0});

    linalg::vector c({5, 1, 1});

    EXPECT_EQ(a + b, c);
}

TEST(VectorOperations, VectorSubtraction) {
    linalg::vector a({2, 0, 1});
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

TEST(VectorMultiplication, IdentityMat) {
    linalg::matrix i2({{1, 0}, {0, 1}});

    linalg::vector a({3, 4});

    EXPECT_EQ(i2 * a, a);
}

TEST(VectorProperties, Magnitude) {
    linalg::vector a({1, 1, 1});
    linalg::vector b({4, 0, 0});
    linalg::vector c({-1, -1, -1});

    EXPECT_DOUBLE_EQ(a.magnitude(), sqrt(3));
    EXPECT_DOUBLE_EQ(b.magnitude(), 4);
    EXPECT_DOUBLE_EQ(a.magnitude(), c.magnitude());
}

TEST(VectorProperties, Normalization) {
    linalg::vector a({7, 3, 2.5});
    linalg::vector a1(a);
    a1.normalize();

    EXPECT_DOUBLE_EQ(a1.magnitude(), 1);
    EXPECT_EQ(a.magnitude() * a1, a);
}

TEST(DotProduct, OrthoDotProd) {
    linalg::vector a({1, 0, 1});
    linalg::vector b({0, 1, 0});

    EXPECT_DOUBLE_EQ(a.dot_product(b).value(), 0);
}

TEST(DotProduct, EqualDotProd) {
    linalg::vector a({2, 1, 1});
    linalg::vector b({2, 1, 1});

    EXPECT_DOUBLE_EQ(a.dot_product(b).value(), a.magnitude() * a.magnitude());
}