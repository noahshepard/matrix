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

} // namespace linalg