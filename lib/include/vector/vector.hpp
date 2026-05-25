#pragma once
#include "matrix/matrix.hpp"
#include <vector>

namespace linalg {
class vector {
public:
  vector(size_t n) : data_(n, 0.0), size_(n) {}
  vector() : data_(), size_(0) {}

  vector(const std::vector<double> &vec) : data_(vec), size_(vec.size()) {}
  vector(std::vector<double> &&vec)
      : data_(std::move(vec)), size_(data_.size()) {}

  vector(const vector &other) : data_(other.data_), size_(other.size_) {}
  vector(vector &&other) : data_(std::move(other.data_)), size_(other.size_) {}

  double magnitude() const;
  double normalize();

  size_t size() const;

  std::optional<double> dot_product(const vector &other) const;

  bool operator==(const vector &other) const;
  bool operator!=(const vector &other) const;

  vector operator+(const vector &other) const;
  vector operator-(const vector &other) const;

  double &operator[](size_t i);
  const double &operator[](size_t i) const;

private:
  friend vector operator*(const vector &v, double f);
  friend vector operator*(double f, const vector &v);

  friend vector operator/(const vector &v, double f);

  friend vector operator*(const matrix &mat, const vector &vec);

  std::vector<double> data_;
  size_t size_;
};

} // namespace linalg