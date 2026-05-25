#include "vector/vector.hpp"
#include "constants.hpp"

namespace linalg {

double vector::magnitude() const {
  double sum = 0;
  for (const double &n : this->data_) {
    sum += n * n;
  }
  return sqrt(sum);
}

double vector::normalize() {
  const double mag = this->magnitude();

  for (double &n : this->data_) {
    n /= mag;
  }

  return mag;
}

size_t vector::size() const { return size_; }

std::optional<double> vector::dot_product(const vector &other) const {
  if (this->size_ != other.size_) {
    return std::nullopt;
  }

  double prod = 0;

  for (size_t i = 0; i < this->size_; ++i) {
    prod += this->data_[i] * other.data_[i];
  }

  return {prod};
}

bool vector::operator==(const vector &other) const {
  if (this->size_ != other.size_) {
    return false;
  }

  for (size_t i = 0; i < this->size_; ++i) {
    if (abs(this->data_[i] - other.data_[i]) > EPS) {
      return false;
    }
  }

  return true;
}

bool vector::operator!=(const vector &other) const { return !(*this == other); }

vector vector::operator+(const vector &other) const {
  size_t sz = (this->size_ < other.size_) ? this->size_ : other.size_;

  vector v(sz);

  for (size_t i = 0; i < sz; i++) {
    v.data_[i] = this->data_[i] + other.data_[i];
  }

  return v;
}

vector vector::operator-(const vector &other) const {
  size_t sz = (this->size_ < other.size_) ? this->size_ : other.size_;

  vector v(sz);

  for (size_t i = 0; i < sz; i++) {
    v.data_[i] = this->data_[i] - other.data_[i];
  }

  return v;
}

double &vector::operator[](size_t i) { return this->data_[i]; }

const double &vector::operator[](size_t i) const { return this->data_[i]; }

vector operator*(const vector &v, double f) {
  vector ret = {};
  ret.data_.reserve(v.size_);
  ret.size_ = v.size_;

  for (const auto &d : v.data_) {
    ret.data_.push_back(d * f);
  }
  return ret;
}

vector operator*(double f, const vector &v) { return operator*(v, f); }

vector operator/(const vector &v, double f) { return operator*(v, 1 / f); }

vector operator*(const matrix &mat, const vector &vec) {

  if (vec.size_ != mat.cols()) {
    throw std::invalid_argument("Wrong transformation size");
  }

  size_t n = vec.size_;
  size_t m = mat.rows();

  vector ret(m);
  // m x n, n x 1
  // n -> m

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      ret.data_[i] += (mat(i, j) * vec.data_[j]);
    }
  }

  return ret;
}

} // namespace linalg