#include "matrix/matrix.hpp"
#include "constants.hpp"
#include <iomanip>

Matrix::Matrix(size_t rows, size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {
  if (rows == 0 || cols == 0)
    throw std::invalid_argument("Matrix dimensions must be >0");
}

Matrix::Matrix(const std::vector<std::vector<double>> &values) {
  if (values.size() == 0 || values[0].size() == 0)
    throw std::invalid_argument("Matrix dimensions must be >0");
  rows_ = values.size();
  cols_ = values[0].size();

  data_.resize(rows_ * cols_);

  for (size_t i = 0; i < rows_; ++i) {
    if (values[i].size() != cols_) {
      throw std::invalid_argument("Matrix must be rectangular");
    }

    for (size_t j = 0; j < cols_; ++j) {
      data_[i * cols_ + j] = values[i][j];
    }
  }
}

Matrix::Matrix(size_t rows, size_t cols, const std::vector<double> &values) {
  if (rows <= 0 || cols <= 0 || values.size() == 0) {
    throw std::invalid_argument("Matrix dimensions must be >0");
  }

  if (rows != values.size() / cols) {
    throw std::invalid_argument("Matrix must be rectangular");
  }

  rows_ = rows;
  cols_ = cols;

  data_.resize(rows_ * cols_);

  for (size_t i = 0; i < rows * cols; ++i) {
    data_[i] = values[i];
  }
}

Matrix::Matrix(size_t n) {
  data_.assign(n * n, 0.0);

  rows_ = n;
  cols_ = n;

  for (size_t r = 0; r < n; r++) {
    data_[r * n + r] = 1.0;
  }
}

double &Matrix::operator()(size_t rows, size_t cols) {
  if (rows >= rows_ || cols >= cols_ || rows < 0 || cols < 0)
    throw std::invalid_argument("Matrix indices out of range");
  return data_[rows * cols_ + cols];
}

double Matrix::operator()(size_t rows, size_t cols) const {
  if (rows >= rows_ || cols >= cols_ || rows < 0 || cols < 0)
    throw std::invalid_argument("Matrix indices out of range");
  return data_[rows * cols_ + cols];
}

bool Matrix::operator==(const Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    return false;

  for (size_t i = 0; i < data_.size(); ++i) {
    if (std::abs(data_[i] - other.data_[i]) > EPS)
      return false;
  }

  return true;
}

bool Matrix::operator!=(const Matrix &other) const { return !(*this == other); }

Matrix Matrix::operator+(const Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Matrix dimensions must match for addition");

  Matrix result(rows_, cols_);

  for (size_t i = 0; i < data_.size(); ++i) {
    result.data_[i] = data_[i] + other.data_[i];
  }

  return result;
}

Matrix Matrix::operator-(const Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Matrix dimensions must match for subtraction");

  Matrix result(rows_, cols_);

  for (size_t i = 0; i < data_.size(); ++i) {
    result.data_[i] = data_[i] - other.data_[i];
  }

  return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
  if (this->cols() != other.rows()) {
    throw std::invalid_argument(
        "The number of cols in the first matrix must match the number of rows "
        "in the second for multiplication");
  }

  Matrix result(this->rows(), other.cols());

  const size_t M = this->rows();
  const size_t P = other.cols();
  const size_t N = this->cols();

  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < P; ++j) {
      double product = 0;
      for (size_t k = 0; k < N; ++k) {
        product += (*this)(i, k) * other(k, j);
      }
      result(i, j) = product;
    }
  }
  return result;
}

size_t Matrix::rows() const { return rows_; }

size_t Matrix::cols() const { return cols_; }

void Matrix::swap_rows(size_t r1, size_t r2) {
  if (r1 >= rows_ || r1 < 0 || r2 >= rows_ || r2 < 0)
    throw std::invalid_argument("Row out of matrix range");

  auto r1_begin = data_.begin() + r1 * cols_;
  auto r2_begin = data_.begin() + r2 * cols_;

  std::swap_ranges(r1_begin, r1_begin + cols_, r2_begin);
}

void Matrix::scale_row(size_t r, double scalar) {
  if (r >= rows_ || r < 0)
    throw std::invalid_argument("Row out of matrix range");

  auto row_begin = data_.begin() + r * cols_;

  for (auto it = row_begin; it != row_begin + cols_; ++it)
    *it *= scalar;
}

void Matrix::add_row_multiple(size_t src, size_t dst, double scalar) {
  if (src >= rows_ || src < 0)
    throw std::invalid_argument("Source out of matrix range");

  if (dst >= rows_ || dst < 0)
    throw std::invalid_argument("Destination out of matrix range");

  auto src_begin = data_.begin() + src * cols_;
  auto dst_begin = data_.begin() + dst * cols_;

  for (size_t i = 0; i < cols_; ++i) {
    dst_begin[i] += scalar * src_begin[i];
  }
}

std::optional<Matrix> Matrix::inverse() const {
  if (this->rows() != this->cols()) {
    return std::nullopt;
  }
  Matrix copy(*this);
  Matrix inv(rows_);

  size_t lead = 0;

  for (size_t row = 0; row < rows_; ++row) {
    size_t pivot_row = row;
    while (pivot_row < rows_ && std::abs((copy)(pivot_row, lead)) < EPS)
      ++pivot_row;
    if (pivot_row == rows_)
      return std::nullopt;

    copy.swap_rows(pivot_row, row);
    inv.swap_rows(pivot_row, row);

    double pivot = (copy)(row, lead);
    copy.scale_row(row, 1.0 / pivot);
    inv.scale_row(row, 1.0 / pivot);

    for (size_t j = 0; j < rows_; ++j) {
      if (j == row)
        continue;
      double factor = (copy)(j, lead);
      if (std::abs(factor) > EPS) {
        copy.add_row_multiple(row, j, -factor);
        inv.add_row_multiple(row, j, -factor);
      }
    }
    ++lead;
  }

  return {inv};
}

double Matrix::determinant() const {
  if (this->rows() != this->cols()) {
    throw std::invalid_argument("Not a square matrix");
  }
  Matrix copy(*this);
  double det = 1;

  size_t lead = 0;

  for (size_t row = 0; row < rows_; ++row) {
    size_t pivot_row = row;
    while (pivot_row < rows_ && std::abs((copy)(pivot_row, lead)) < EPS)
      ++pivot_row;
    if (pivot_row == rows_)
      return 0;

    if (pivot_row != row) {
      copy.swap_rows(pivot_row, row);
      det *= -1;
    }

    double pivot = (copy)(row, lead);
    copy.scale_row(row, 1.0 / pivot);
    det *= (pivot);

    for (size_t j = pivot_row + 1; j < rows_; ++j) {
      double factor = (copy)(j, lead);
      if (std::abs(factor) > EPS) {
        copy.add_row_multiple(row, j, -factor);
      }
    }
    ++lead;
  }

  return det;
}

void Matrix::rref() {
  size_t lead = 0;

  for (size_t row = 0; row < rows_; ++row) {
    if (lead >= cols_)
      return;

    size_t i = row;
    while (i < rows_ && std::abs((*this)(i, lead)) < EPS) {
      ++i;
    }

    if (i == rows_) {
      lead++;
      row--;
      continue;
    }

    swap_rows(i, row);

    double pivot = (*this)(row, lead);
    if (std::abs(pivot) > EPS) {
      scale_row(row, 1.0 / pivot);
    }

    for (size_t j = 0; j < rows_; ++j) {
      if (j != row) {
        double factor = (*this)(j, lead);
        if (std::abs(factor) > EPS)
          add_row_multiple(row, j, -factor);
      }
    }

    ++lead;
  }
}

std::ostream &operator<<(std::ostream &os, const Matrix &m) {
  for (size_t r = 0; r < m.rows_; ++r) {
    os << "[ ";
    for (size_t c = 0; c < m.cols_; ++c) {
      os << std::setw(8) << m.data_[r * m.cols_ + c] << " ";
    }
    os << "]\n";
  }
  return os;
}

Matrix operator*(const Matrix &m, double f) {
  Matrix result(m.rows_, m.cols_);
  for (size_t i = 0; i < m.data_.size(); ++i) {
    result.data_[i] = m.data_[i] * f;
  }
  return result;
}

Matrix operator*(double f, const Matrix &m) { return m * f; }

Matrix operator/(const Matrix &m, double f) {
  if (std::abs(f) < EPS) {
    throw std::invalid_argument("Division by zero");
  }

  return m * (1.0 / f);
}
