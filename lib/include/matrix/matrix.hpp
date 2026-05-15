#pragma once

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

namespace linalg {

class matrix {
public:
  matrix(size_t rows, size_t cols);

  matrix(const std::vector<std::vector<double>> &values);

  matrix(size_t rows, size_t cols, const std::vector<double> &values);

  matrix(const matrix &other) noexcept
      : data_(other.data_), rows_(other.rows_), cols_(other.cols_) {}

  matrix(matrix &&other) noexcept
      : data_(std::move(other.data_)), rows_(other.rows_), cols_(other.cols_) {}

  static matrix identity(size_t n);

  double &operator()(size_t rows, size_t cols);
  double operator()(size_t rows, size_t cols) const;

  bool operator==(const matrix &other) const;
  bool operator!=(const matrix &other) const;

  matrix operator+(const matrix &other) const;
  matrix operator-(const matrix &other) const;

  matrix operator*(const matrix &other) const;

  size_t rows() const;
  size_t cols() const;

  std::optional<matrix> inverse() const;
  double determinant() const;

  void rref();

private:
  void swap_rows(size_t r1, size_t r2);
  void scale_row(size_t r, double scalar);
  void add_row_multiple(size_t src, size_t dst, double scalar);

  friend std::ostream &operator<<(std::ostream &os, const matrix &m);

  friend matrix operator*(const matrix &m, double f);
  friend matrix operator*(double f, const matrix &m);

  friend matrix operator/(const matrix &m, double f);

  size_t rows_;
  size_t cols_;
  std::vector<double> data_;
};

matrix identity(size_t n);

} // namespace linalg
