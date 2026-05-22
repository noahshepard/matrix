#include "vector.hpp"

namespace linalg {

vector::vector(const std::vector<int> &vec) {
  this->size_ = vec.size();
  this->data_.reserve(this->size_);
  for (const int &num : vec) {
    this->data_.push_back(static_cast<double>(num));
  }
}

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

  for (const auto &d : v.data_) {
    ret.data_.push_back(d * f);
  }
  return ret;
}

vector operator*(double f, const vector &v) { return operator*(v, f); }

vector operator/(const vector &v, double f) { return operator*(v, 1 / f); }

} // namespace linalg