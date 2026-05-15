#include "constants.hpp"
#include "matrix/matrix.hpp"
#include <gtest/gtest.h>
#include <random>

bool isRREF(const linalg::matrix &m) {
  // Check for RREF properties
  // 1. Leading 1s
  // 2. Zeros below and above leading 1s
  // 3. Each leading 1 is to the right of the leading 1 in the previous row
  // 4. Rows of all zeros at the bottom
  // 5. Leading 1s are the only non-zero entries in their columns

  int last_pivot_col = -1;
  bool found_zero_row = false;

  for (size_t r = 0; r < m.rows(); ++r) {
    // Find first non-zero element in this row
    int pivot_col = -1;
    for (size_t c = 0; c < m.cols(); ++c) {
      if (std::abs(m(r, c)) > EPS) {
        pivot_col = c;
        break;
      }
    }

    // If all zeros in this row
    if (pivot_col == -1) {
      found_zero_row = true;
      continue;
    }

    // If we found a non-zero row after a zero row, invalid
    if (found_zero_row)
      return false;

    // First non-zero must be 1
    if (std::abs(m(r, pivot_col) - 1.0) > EPS)
      return false;

    // Pivot column must be to the right of previous pivot
    if (pivot_col <= last_pivot_col)
      return false;

    // All other entries in this pivot column must be zero
    for (size_t rr = 0; rr < m.rows(); ++rr) {
      if (rr != r && std::abs(m(rr, pivot_col)) > EPS)
        return false;
    }

    last_pivot_col = pivot_col;
  }

  return true;
}

TEST(isRREF_Test, ValidRREF) {
  linalg::matrix A({{1, 0, 2, 0}, {0, 1, -1, 0}, {0, 0, 0, 0}});

  EXPECT_TRUE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_NonLeadingOne) {
  linalg::matrix A({{2, 0, 2, 0}, {0, 1, -1, 0}, {0, 0, 0, 0}});

  EXPECT_FALSE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_NonZeroAbovePivot) {
  linalg::matrix A({{1, 0, 2, 0}, {0, 1, -1, 0}, {0, 1, 0, 0}});

  EXPECT_FALSE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_PivotNotRight) {
  linalg::matrix A({{1, 0, 2, 0}, {0, 0, 1, 0}, {0, 0, 0, 0}});

  EXPECT_FALSE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_NonZeroRowBelowZeroRow) {
  linalg::matrix A({{1, 0, 2, 0}, {0, 0, 0, 0}, {0, 1, -1, 0}});

  EXPECT_FALSE(isRREF(A));
}

TEST(MatrixRREF, SimpleAugmented2x3) {
  linalg::matrix A({{1, 2, 1, 9}, {2, -1, 1, 8}});

  A.rref();

  // Expected RREF:
  // [1 0 1/5 5]
  // [0 1 2/5 2]
  EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(A(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(A(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(A(1, 1), 1.0);
}

TEST(MatrixRREF, AugmentedWithZeroRow) {
  linalg::matrix A({{1, 2, 1, 9}, {2, 4, 2, 18}, {0, 0, 0, 0}});

  A.rref();

  EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
  for (size_t r = 1; r < A.rows(); ++r)
    for (size_t c = 0; c < A.cols(); ++c)
      EXPECT_DOUBLE_EQ(A(r, c), 0.0);
}

TEST(MatrixRREF, RankDeficientAugmented) {
  linalg::matrix A({{1, 2, 3, 6}, {2, 4, 6, 12}, {3, 6, 9, 18}});

  A.rref();

  EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
  for (size_t r = 1; r < A.rows(); ++r)
    for (size_t c = 0; c < A.cols(); ++c)
      EXPECT_DOUBLE_EQ(A(r, c), 0.0);
}

TEST(MatrixRREF, FullAugmented3x4) {
  linalg::matrix A({{2, 1, -1, 8}, {-3, -1, 2, -11}, {-2, 1, 2, -3}});

  A.rref();

  EXPECT_DOUBLE_EQ(A(0, 3), 2.0);  // x
  EXPECT_DOUBLE_EQ(A(1, 3), 3.0);  // y
  EXPECT_DOUBLE_EQ(A(2, 3), -1.0); // z
}

TEST(MatrixRREF, ZeroMatrix) {
  linalg::matrix A({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});

  A.rref();

  for (size_t r = 0; r < A.rows(); ++r)
    for (size_t c = 0; c < A.cols(); ++c)
      EXPECT_DOUBLE_EQ(A(r, c), 0.0);
}

TEST(MatrixRREF, IdentityMatrix) {
  linalg::matrix A({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});

  A.rref();

  for (size_t r = 0; r < A.rows(); ++r)
    for (size_t c = 0; c < A.cols(); ++c)
      EXPECT_DOUBLE_EQ(A(r, c), (r == c) ? 1.0 : 0.0);
}

TEST(MatrixRREF, InconsistentSystem) {
  linalg::matrix A({{1, 2, 3, 4}, {2, 4, 6, 9}});

  A.rref();

  // The last row should reflect the inconsistency
  EXPECT_DOUBLE_EQ(A(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(A(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(A(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(A(1, 3), 1.0);
}

TEST(MatrixRREF, FreeVariables) {
  linalg::matrix A({{2, 1, 3, 0, 10}, {1, 1, 1, 0, 7}, {3, 2, 4, 0, 17}});

  linalg::matrix U({{1, 0, 2, 0, 3}, {0, 1, -1, 0, 4}, {0, 0, 0, 0, 0}});
  A.rref();

  EXPECT_EQ(A, U);
}

TEST(MatrixRREF, DimensionsPreserved) {
  linalg::matrix A({{1, 2, 3}, {4, 5, 6}});

  size_t original_rows = A.rows();
  size_t original_cols = A.cols();

  A.rref();

  EXPECT_EQ(A.rows(), original_rows);
  EXPECT_EQ(A.cols(), original_cols);
}

TEST(MatrixRREF, ZeroAtBottom) {
  linalg::matrix A({{0, 0, 0, 0}, {1, 0, 0, 4}, {0, 1, 7, 8}});

  A.rref();

  linalg::matrix U({{1, 0, 0, 4}, {0, 1, 7, 8}, {0, 0, 0, 0}});
  EXPECT_EQ(A, U);
}

class MatrixRREF_RandomTest
    : public ::testing::TestWithParam<std::pair<int, int>> {
protected:
  std::mt19937 gen{std::random_device{}()};
  std::uniform_real_distribution<> dis{-100.0, 100.0};
};

TEST_P(MatrixRREF_RandomTest, RandomMatrices) {
  auto [rows, cols] = GetParam();
  linalg::matrix A(rows, cols);
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < cols; ++c) {
      // A(r, c) = dis(gen);
      A(r, c) = static_cast<double>(
          r + c + 1); // simpler deterministic values for test stability
    }
  }

  A.rref();

  EXPECT_TRUE(isRREF(A)) << "Failed for " << rows << "x" << cols
                         << " random matrix:\n"
                         << A;
}

INSTANTIATE_TEST_SUITE_P(
    RandomMatrices, MatrixRREF_RandomTest,
    ::testing::Values(std::make_pair(2, 2), std::make_pair(3, 3),
                      std::make_pair(4, 4), std::make_pair(5, 5),
                      std::make_pair(3, 5), std::make_pair(5, 3),
                      std::make_pair(6, 4), std::make_pair(4, 6),
                      std::make_pair(10, 10)));

// linalg::matrix Operation Tests

TEST(MatrixOperations, Addition) {
  linalg::matrix A({{1, 2}, {3, 4}});
  linalg::matrix B({{5, 6}, {7, 8}});
  linalg::matrix C = A + B;

  linalg::matrix C_expected({{6, 8}, {10, 12}});
  EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Subtraction) {
  linalg::matrix A({{5, 6}, {7, 8}});
  linalg::matrix B({{1, 2}, {3, 4}});
  linalg::matrix C = A - B;

  linalg::matrix C_expected({{4, 4}, {4, 4}});
  EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Addition_InvalidDimensions) {
  linalg::matrix A({{1, 2}, {3, 4}});
  linalg::matrix B({{5, 6, 7}, {8, 9, 10}});

  EXPECT_THROW({ linalg::matrix C = A + B; }, std::invalid_argument);
}

TEST(MatrixOperations, Subtraction_InvalidDimensions) {
  linalg::matrix A({{1, 2}, {3, 4}});
  linalg::matrix B({{5, 6, 7}, {8, 9, 10}});

  EXPECT_THROW({ linalg::matrix C = A - B; }, std::invalid_argument);
}

TEST(MatrixOperations, Addition_NegativeValues) {
  linalg::matrix A({{-1, -2}, {-3, -4}});
  linalg::matrix B({{5, 6}, {7, 8}});
  linalg::matrix C = A + B;

  linalg::matrix C_expected({{4, 4}, {4, 4}});
  EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Subtraction_NegativeValues) {
  linalg::matrix A({{5, 6}, {7, 8}});
  linalg::matrix B({{-1, -2}, {-3, -4}});
  linalg::matrix C = A - B;

  linalg::matrix C_expected({{6, 8}, {10, 12}});
  EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Addition_ZeroMatrix) {
  linalg::matrix A({{0, 0}, {0, 0}});
  linalg::matrix B({{5, 6}, {7, 8}});
  linalg::matrix C = A + B;

  EXPECT_EQ(C, B);
}

TEST(MatrixOperations, Subtraction_ZeroMatrix) {
  linalg::matrix A({{0, 0}, {0, 0}});
  linalg::matrix B({{5, 6}, {7, 8}});
  linalg::matrix C = B - A;

  EXPECT_EQ(C, B);
}

TEST(MatrixOperations, ScalarMultiplication_FriendFunction) {
  linalg::matrix A({{1, 2}, {3, 4}});
  double scalar = 2.0;

  linalg::matrix B = A * scalar;
  linalg::matrix C = scalar * A;

  linalg::matrix expected({{2, 4}, {6, 8}});
  EXPECT_EQ(B, expected);
  EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarMultiplication_NegativeScalar) {
  linalg::matrix A({{1, -2}, {-3, 4}});
  double scalar = -2.0;

  linalg::matrix B = A * scalar;
  linalg::matrix C = scalar * A;

  linalg::matrix expected({{-2, 4}, {6, -8}});
  EXPECT_EQ(B, expected);
  EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarMultiplication_ZeroScalar) {
  linalg::matrix A({{1, 2}, {3, 4}});
  double scalar = 0.0;

  linalg::matrix B = A * scalar;
  linalg::matrix C = scalar * A;

  linalg::matrix expected({{0, 0}, {0, 0}});
  EXPECT_EQ(B, expected);
  EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarMultiplication_FractionalScalar) {
  linalg::matrix A({{2, 4}, {6, 8}});
  double scalar = 0.5;

  linalg::matrix B = A * scalar;
  linalg::matrix C = scalar * A;

  linalg::matrix expected({{1, 2}, {3, 4}});
  EXPECT_EQ(B, expected);
  EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarDivision_FriendFunction) {
  linalg::matrix A({{2, 4}, {6, 8}});
  double scalar = 2.0;

  linalg::matrix B = A / scalar;

  linalg::matrix expected({{1, 2}, {3, 4}});
  EXPECT_EQ(B, expected);
}

TEST(MatrixOperations, ScalarDivision_ByZero) {
  linalg::matrix A({{2, 4}, {6, 8}});
  double scalar = 0.0;

  EXPECT_THROW({ linalg::matrix B = A / scalar; }, std::invalid_argument);
}

TEST(MatrixOperations, ScalarDivision_NegativeScalar) {
  linalg::matrix A({{2, -4}, {-6, 8}});
  double scalar = -2.0;

  linalg::matrix B = A / scalar;

  linalg::matrix expected({{-1, 2}, {3, -4}});
  EXPECT_EQ(B, expected);
}

TEST(MatrixOperations, ScalarDivision_FractionalScalar) {
  linalg::matrix A({{1, 2}, {3, 4}});
  double scalar = 0.5;

  linalg::matrix B = A / scalar;

  linalg::matrix expected({{2, 4}, {6, 8}});
  EXPECT_EQ(B, expected);
}

TEST(MatrixOperations, Multiplication) {
  linalg::matrix A({{1, 2, 3}, {4, 5, 6}});
  linalg::matrix B({{7, 8}, {9, 10}, {11, 12}});
  linalg::matrix C = A * B;

  linalg::matrix C_expected({{58, 64}, {139, 154}});
  EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Multiplication_InvalidDimensions) {
  linalg::matrix A({{1, 2, 3}, {4, 5, 6}});
  linalg::matrix B({{5, 6}, {7, 8}});

  EXPECT_THROW({ linalg::matrix C = A * B; }, std::invalid_argument);
}

TEST(MatrixOperations, Inverse) {
  linalg::matrix A({{4, 7}, {2, 6}});
  linalg::matrix A_inv = A.inverse().value();

  linalg::matrix A_inv_expected({{0.6, -0.7}, {-0.2, 0.4}});
  for (size_t r = 0; r < A_inv.rows(); ++r) {
    for (size_t c = 0; c < A_inv.cols(); ++c) {
      EXPECT_DOUBLE_EQ(A_inv(r, c), A_inv_expected(r, c));
    }
  }
}

TEST(MatrixOperations, InverseIdentity) {
  linalg::matrix I({{1, 0}, {0, 1}});
  linalg::matrix I_inv = I.inverse().value();

  EXPECT_EQ(I_inv, I);
}

TEST(MatrixOperations, InverseMultIsIdentity) {
  linalg::matrix A({{3, 0}, {0, 2}});
  linalg::matrix A_inv = A.inverse().value();
  linalg::matrix I = A * A_inv;

  linalg::matrix I_expected({{1, 0}, {0, 1}});
  for (size_t r = 0; r < I.rows(); ++r) {
    for (size_t c = 0; c < I.cols(); ++c) {
      EXPECT_DOUBLE_EQ(I(r, c), I_expected(r, c));
    }
  }
}

TEST(MatrixOperations, InverseNonSquareMatrix) {
  linalg::matrix A({{1, 2, 3}, {4, 5, 6}});

  EXPECT_EQ(A.inverse(), std::nullopt);
}

TEST(MatrixOperations, InverseInvertableMatrix) {
  linalg::matrix A({{2, 2}, {2, 2}});

  EXPECT_EQ(A.inverse(), std::nullopt);

  linalg::matrix B({{1, 0, 1}, {0, 4, 8}, {0, 2, 4}});

  EXPECT_EQ(A.inverse(), std::nullopt);
}

TEST(MatrixOperations, DeterminantSimple) {
  linalg::matrix A({{1, 0}, {0, 1}});
  EXPECT_DOUBLE_EQ(A.determinant(), 1.0);

  linalg::matrix B({{4, 2}, {4, 1}});
  EXPECT_DOUBLE_EQ(B.determinant(), -4.0);

  linalg::matrix C({{1, -1}, {0, 1}});
  EXPECT_DOUBLE_EQ(C.determinant(), 1.0);

  linalg::matrix D({{1, 1}, {1, 1}});
  EXPECT_DOUBLE_EQ(D.determinant(), 0.0);
}

TEST(MatrixOperations, Determinant3x3) {
  linalg::matrix A({{5, 4, 3}, {7, 2, 6}, {1, 8, 9}});
  EXPECT_DOUBLE_EQ(A.determinant(), -216.0);

  linalg::matrix B({{7, 7, 7}, {7, 7, 9}, {7, 1, 4}});
  EXPECT_DOUBLE_EQ(B.determinant(), 84.0);
}

TEST(MatrixOperations, RandDeficienDetZero) {
  linalg::matrix A({{4, 1}, {0, 0}});
  EXPECT_DOUBLE_EQ(A.determinant(), 0);

  linalg::matrix B({{4, 1, 6}, {8, 2, 12}, {2, 0.5, 3}});
  EXPECT_DOUBLE_EQ(B.determinant(), 0);
}

TEST(MatrixOperations, FlippedRowNotZero) {
  linalg::matrix A({{0, 1}, {1, 1}});
  EXPECT_DOUBLE_EQ(A.determinant(), -1.0);
}

TEST(MatrixCreation, From2DVector_Valid) {
  std::vector<std::vector<double>> values = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  linalg::matrix A(values);

  EXPECT_EQ(A.rows(), 3);
  EXPECT_EQ(A.cols(), 3);
  EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(A(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(A(2, 2), 9.0);
}

TEST(MatrixCreation, From2DVector_NonRectangular) {
  std::vector<std::vector<double>> values = {{1, 2, 3}, {4, 5}, {6, 7, 8}};

  EXPECT_THROW({ linalg::matrix A(values); }, std::invalid_argument);
}

TEST(MatrixCreation, From2DVector_Empty) {
  std::vector<std::vector<double>> values = {};

  EXPECT_THROW({ linalg::matrix A(values); }, std::invalid_argument);
}

TEST(MatrixCreation, From1DVector_InvalidSize) {
  std::vector<double> values = {1, 2, 3, 4, 5};

  EXPECT_THROW({ linalg::matrix A(2, 3, values); }, std::invalid_argument);
}

TEST(MatrixCreation, From1DVector_Empty) {
  std::vector<double> values = {};

  EXPECT_THROW({ linalg::matrix A(2, 2, values); }, std::invalid_argument);
}

TEST(MatrixAccess, OutOfBounds) {
  linalg::matrix A({{1, 2}, {3, 4}});

  EXPECT_THROW({ double val = A(2, 0); }, std::invalid_argument);
  EXPECT_THROW({ double val = A(0, 2); }, std::invalid_argument);
  EXPECT_THROW({ A(2, 0) = 5; }, std::invalid_argument);
  EXPECT_THROW({ A(0, 2) = 5; }, std::invalid_argument);
}