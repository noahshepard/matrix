#include <gtest/gtest.h>
#include "matrix/matrix.hpp"
#include "constants.hpp"
#include <random>

// Helper to compare doubles with EPS tolerance
void EXPECT_DOUBLE_EQ_EPS(double a, double b)
{
    EXPECT_NEAR(a, b, EPS);
}

bool isRREF(const Matrix &m)
{
    // Check for RREF properties
    // 1. Leading 1s
    // 2. Zeros below and above leading 1s
    // 3. Each leading 1 is to the right of the leading 1 in the previous row
    // 4. Rows of all zeros at the bottom
    // 5. Leading 1s are the only non-zero entries in their columns

    int last_pivot_col = -1;
    bool found_zero_row = false;

    for (size_t r = 0; r < m.rows(); ++r)
    {
        // Find first non-zero element in this row
        int pivot_col = -1;
        for (size_t c = 0; c < m.cols(); ++c)
        {
            if (std::abs(m(r, c)) > EPS)
            {
                pivot_col = c;
                break;
            }
        }

        // If all zeros in this row
        if (pivot_col == -1)
        {
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
        for (size_t rr = 0; rr < m.rows(); ++rr)
        {
            if (rr != r && std::abs(m(rr, pivot_col)) > EPS)
                return false;
        }

        last_pivot_col = pivot_col;
    }

    return true;
}

TEST(isRREF_Test, ValidRREF)
{
    Matrix A({{1, 0, 2, 0},
              {0, 1, -1, 0},
              {0, 0, 0, 0}});

    EXPECT_TRUE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_NonLeadingOne)
{
    Matrix A({{2, 0, 2, 0},
              {0, 1, -1, 0},
              {0, 0, 0, 0}});

    EXPECT_FALSE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_NonZeroAbovePivot)
{
    Matrix A({{1, 0, 2, 0},
              {0, 1, -1, 0},
              {0, 1, 0, 0}});

    EXPECT_FALSE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_PivotNotRight)
{
    Matrix A({{1, 0, 2, 0},
              {0, 0, 1, 0},
              {0, 0, 0, 0}});

    EXPECT_FALSE(isRREF(A));
}

TEST(isRREF_Test, InvalidRREF_NonZeroRowBelowZeroRow)
{
    Matrix A({{1, 0, 2, 0},
              {0, 0, 0, 0},
              {0, 1, -1, 0}});

    EXPECT_FALSE(isRREF(A));
}

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

TEST(MatrixRREF, AugmentedWithZeroRow)
{
    Matrix A({{1, 2, 1, 9},
              {2, 4, 2, 18},
              {0, 0, 0, 0}});

    A.rref();

    EXPECT_DOUBLE_EQ_EPS(A(0, 0), 1.0);
    for (size_t r = 1; r < A.rows(); ++r)
        for (size_t c = 0; c < A.cols(); ++c)
            EXPECT_DOUBLE_EQ_EPS(A(r, c), 0.0);
}

TEST(MatrixRREF, RankDeficientAugmented)
{
    Matrix A({{1, 2, 3, 6},
              {2, 4, 6, 12},
              {3, 6, 9, 18}});

    A.rref();

    EXPECT_DOUBLE_EQ_EPS(A(0, 0), 1.0);
    for (size_t r = 1; r < A.rows(); ++r)
        for (size_t c = 0; c < A.cols(); ++c)
            EXPECT_DOUBLE_EQ_EPS(A(r, c), 0.0);
}

TEST(MatrixRREF, FullAugmented3x4)
{
    Matrix A({{2, 1, -1, 8},
              {-3, -1, 2, -11},
              {-2, 1, 2, -3}});

    A.rref();

    EXPECT_DOUBLE_EQ_EPS(A(0, 3), 2.0);  // x
    EXPECT_DOUBLE_EQ_EPS(A(1, 3), 3.0);  // y
    EXPECT_DOUBLE_EQ_EPS(A(2, 3), -1.0); // z
}

TEST(MatrixRREF, ZeroMatrix)
{
    Matrix A({{0, 0, 0},
              {0, 0, 0},
              {0, 0, 0}});

    A.rref();

    for (size_t r = 0; r < A.rows(); ++r)
        for (size_t c = 0; c < A.cols(); ++c)
            EXPECT_DOUBLE_EQ_EPS(A(r, c), 0.0);
}

TEST(MatrixRREF, IdentityMatrix)
{
    Matrix A({{1, 0, 0},
              {0, 1, 0},
              {0, 0, 1}});

    A.rref();

    for (size_t r = 0; r < A.rows(); ++r)
        for (size_t c = 0; c < A.cols(); ++c)
            EXPECT_DOUBLE_EQ_EPS(A(r, c), (r == c) ? 1.0 : 0.0);
}

TEST(MatrixRREF, InconsistentSystem)
{
    Matrix A({{1, 2, 3, 4},
              {2, 4, 6, 9}});

    A.rref();

    // The last row should reflect the inconsistency
    EXPECT_DOUBLE_EQ_EPS(A(1, 0), 0.0);
    EXPECT_DOUBLE_EQ_EPS(A(1, 1), 0.0);
    EXPECT_DOUBLE_EQ_EPS(A(1, 2), 0.0);
    EXPECT_DOUBLE_EQ_EPS(A(1, 3), 1.0);
}

TEST(MatrixRREF, FreeVariables)
{
    Matrix A({{2, 1, 3, 0, 10},
              {1, 1, 1, 0, 7},
              {3, 2, 4, 0, 17}});

    Matrix U({{1, 0, 2, 0, 3},
              {0, 1, -1, 0, 4},
              {0, 0, 0, 0, 0}});
    A.rref();

    EXPECT_EQ(A, U);
}

TEST(MatrixRREF, DimensionsPreserved)
{
    Matrix A({{1, 2, 3},
              {4, 5, 6}});

    size_t original_rows = A.rows();
    size_t original_cols = A.cols();

    A.rref();

    EXPECT_EQ(A.rows(), original_rows);
    EXPECT_EQ(A.cols(), original_cols);
}

TEST(MatrixRREF, ZeroAtBottom)
{
    Matrix A({{0, 0, 0, 0},
              {1, 0, 0, 4},
              {0, 1, 7, 8}});

    A.rref();

    Matrix U({{1, 0, 0, 4},
              {0, 1, 7, 8},
              {0, 0, 0, 0}});
    EXPECT_EQ(A, U);
}

class MatrixRREF_RandomTest : public ::testing::TestWithParam<std::pair<int, int>>
{
protected:
    std::mt19937 gen{std::random_device{}()};
    std::uniform_real_distribution<> dis{-100.0, 100.0};
};

TEST_P(MatrixRREF_RandomTest, RandomMatrices)
{
    auto [rows, cols] = GetParam();
    Matrix A(rows, cols);
    for (size_t r = 0; r < rows; ++r)
    {
        for (size_t c = 0; c < cols; ++c)
        {
            // A(r, c) = dis(gen);
            A(r, c) = static_cast<double>(r + c + 1); // simpler deterministic values for test stability
        }
    }

    A.rref();

    EXPECT_TRUE(isRREF(A)) << "Failed for " << rows << "x" << cols << " random matrix:\n"
                           << A;
}

INSTANTIATE_TEST_SUITE_P(
    RandomMatrices,
    MatrixRREF_RandomTest,
    ::testing::Values(
        std::make_pair(2, 2),
        std::make_pair(3, 3),
        std::make_pair(4, 4),
        std::make_pair(5, 5),
        std::make_pair(3, 5),
        std::make_pair(5, 3),
        std::make_pair(6, 4),
        std::make_pair(4, 6),
        std::make_pair(10, 10)));

// Matrix Operation Tests

TEST(MatrixOperations, Addition)
{
    Matrix A({{1, 2}, {3, 4}});
    Matrix B({{5, 6}, {7, 8}});
    Matrix C = A + B;

    Matrix C_expected({{6, 8}, {10, 12}});
    EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Subtraction)
{
    Matrix A({{5, 6}, {7, 8}});
    Matrix B({{1, 2}, {3, 4}});
    Matrix C = A - B;

    Matrix C_expected({{4, 4}, {4, 4}});
    EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Addition_InvalidDimensions)
{
    Matrix A({{1, 2}, {3, 4}});
    Matrix B({{5, 6, 7}, {8, 9, 10}});

    EXPECT_THROW({ Matrix C = A + B; }, std::invalid_argument);
}

TEST(MatrixOperations, Subtraction_InvalidDimensions)
{
    Matrix A({{1, 2}, {3, 4}});
    Matrix B({{5, 6, 7}, {8, 9, 10}});

    EXPECT_THROW({ Matrix C = A - B; }, std::invalid_argument);
}

TEST(MatrixOperations, Addition_NegativeValues)
{
    Matrix A({{-1, -2}, {-3, -4}});
    Matrix B({{5, 6}, {7, 8}});
    Matrix C = A + B;

    Matrix C_expected({{4, 4}, {4, 4}});
    EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Subtraction_NegativeValues)
{
    Matrix A({{5, 6}, {7, 8}});
    Matrix B({{-1, -2}, {-3, -4}});
    Matrix C = A - B;

    Matrix C_expected({{6, 8}, {10, 12}});
    EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Addition_ZeroMatrix)
{
    Matrix A({{0, 0}, {0, 0}});
    Matrix B({{5, 6}, {7, 8}});
    Matrix C = A + B;

    EXPECT_EQ(C, B);
}

TEST(MatrixOperations, Subtraction_ZeroMatrix)
{
    Matrix A({{0, 0}, {0, 0}});
    Matrix B({{5, 6}, {7, 8}});
    Matrix C = B - A;

    EXPECT_EQ(C, B);
}

TEST(MatrixOperations, ScalarMultiplication_FriendFunction)
{
    Matrix A({{1, 2}, {3, 4}});
    double scalar = 2.0;

    Matrix B = A * scalar;
    Matrix C = scalar * A;

    Matrix expected({{2, 4}, {6, 8}});
    EXPECT_EQ(B, expected);
    EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarMultiplication_NegativeScalar)
{
    Matrix A({{1, -2}, {-3, 4}});
    double scalar = -2.0;

    Matrix B = A * scalar;
    Matrix C = scalar * A;

    Matrix expected({{-2, 4}, {6, -8}});
    EXPECT_EQ(B, expected);
    EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarMultiplication_ZeroScalar)
{
    Matrix A({{1, 2}, {3, 4}});
    double scalar = 0.0;

    Matrix B = A * scalar;
    Matrix C = scalar * A;

    Matrix expected({{0, 0}, {0, 0}});
    EXPECT_EQ(B, expected);
    EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarMultiplication_FractionalScalar)
{
    Matrix A({{2, 4}, {6, 8}});
    double scalar = 0.5;

    Matrix B = A * scalar;
    Matrix C = scalar * A;

    Matrix expected({{1, 2}, {3, 4}});
    EXPECT_EQ(B, expected);
    EXPECT_EQ(C, expected);
}

TEST(MatrixOperations, ScalarDivision_FriendFunction)
{
    Matrix A({{2, 4}, {6, 8}});
    double scalar = 2.0;

    Matrix B = A / scalar;

    Matrix expected({{1, 2}, {3, 4}});
    EXPECT_EQ(B, expected);
}

TEST(MatrixOperations, ScalarDivision_ByZero)
{
    Matrix A({{2, 4}, {6, 8}});
    double scalar = 0.0;

    EXPECT_THROW({ Matrix B = A / scalar; }, std::invalid_argument);
}

TEST(MatrixOperations, ScalarDivision_NegativeScalar)
{
    Matrix A({{2, -4}, {-6, 8}});
    double scalar = -2.0;

    Matrix B = A / scalar;

    Matrix expected({{-1, 2}, {3, -4}});
    EXPECT_EQ(B, expected);
}

TEST(MatrixOperations, ScalarDivision_FractionalScalar)
{
    Matrix A({{1, 2}, {3, 4}});
    double scalar = 0.5;

    Matrix B = A / scalar;

    Matrix expected({{2, 4}, {6, 8}});
    EXPECT_EQ(B, expected);
}

TEST(MatrixOperations, Multiplication)
{
    Matrix A({{1, 2, 3},
              {4, 5, 6}});
    Matrix B({{7, 8},
              {9, 10},
              {11, 12}});
    Matrix C = A * B;

    Matrix C_expected({{58, 64},
                       {139, 154}});
    EXPECT_EQ(C, C_expected);
}

TEST(MatrixOperations, Multiplication_InvalidDimensions)
{
    Matrix A({{1, 2, 3}, {4, 5, 6}});
    Matrix B({{5, 6}, {7, 8}});

    EXPECT_THROW({ Matrix C = A * B; }, std::invalid_argument);
}

TEST(MatrixOperations, Inverse)
{
    Matrix A({{4, 7},
              {2, 6}});
    Matrix A_inv = A.inverse();

    Matrix A_inv_expected({{0.6, -0.7},
                           {-0.2, 0.4}});
    for (size_t r = 0; r < A_inv.rows(); ++r)
    {
        for (size_t c = 0; c < A_inv.cols(); ++c)
        {
            EXPECT_DOUBLE_EQ_EPS(A_inv(r, c), A_inv_expected(r, c));
        }
    }
}

TEST(MatrixOperations, InverseIdentity)
{
    Matrix I({{1, 0},
              {0, 1}});
    Matrix I_inv = I.inverse();

    EXPECT_EQ(I_inv, I);
}

TEST(MatrixOperations, InverseMultIsIdentity)
{
    Matrix A({{3, 0},
              {0, 2}});
    Matrix A_inv = A.inverse();
    Matrix I = A * A_inv;

    Matrix I_expected({{1, 0},
                       {0, 1}});
    for (size_t r = 0; r < I.rows(); ++r)
    {
        for (size_t c = 0; c < I.cols(); ++c)
        {
            EXPECT_DOUBLE_EQ_EPS(I(r, c), I_expected(r, c));
        }
    }
}

TEST(MatrixOperations, InverseNonSquareMatrix)
{
    Matrix A({{1, 2, 3},
              {4, 5, 6}});

    EXPECT_THROW({ Matrix A_inv = A.inverse(); }, std::invalid_argument);
}

TEST(MatrixCreation, From2DVector_Valid)
{
    std::vector<std::vector<double>> values = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}};

    Matrix A(values);

    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_DOUBLE_EQ_EPS(A(0, 0), 1.0);
    EXPECT_DOUBLE_EQ_EPS(A(1, 1), 5.0);
    EXPECT_DOUBLE_EQ_EPS(A(2, 2), 9.0);
}

TEST(MatrixCreation, From2DVector_NonRectangular)
{
    std::vector<std::vector<double>> values = {
        {1, 2, 3},
        {4, 5},
        {6, 7, 8}};

    EXPECT_THROW({ Matrix A(values); }, std::invalid_argument);
}

TEST(MatrixCreation, From2DVector_Empty)
{
    std::vector<std::vector<double>> values = {};

    EXPECT_THROW({ Matrix A(values); }, std::invalid_argument);
}

TEST(MatrixCreation, From1DVector_InvalidSize)
{
    std::vector<double> values = {1, 2, 3, 4, 5};

    EXPECT_THROW({ Matrix A(2, 3, values); }, std::invalid_argument);
}

TEST(MatrixCreation, From1DVector_Empty)
{
    std::vector<double> values = {};

    EXPECT_THROW({ Matrix A(2, 2, values); }, std::invalid_argument);
}

TEST(MatrixAccess, OutOfBounds)
{
    Matrix A({{1, 2}, {3, 4}});

    EXPECT_THROW({ double val = A(2, 0); }, std::invalid_argument);
    EXPECT_THROW({ double val = A(0, 2); }, std::invalid_argument);
    EXPECT_THROW({ A(2, 0) = 5; }, std::invalid_argument);
    EXPECT_THROW({ A(0, 2) = 5; }, std::invalid_argument);
}