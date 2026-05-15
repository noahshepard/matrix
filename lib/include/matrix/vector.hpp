#pragma once
#include <vector>

namespace linalg {
class vector {
public:
  vector(size_t n);

  vector(const std::vector<double> &vec);
  vector(const std::vector<int> &vec);

  vector(std::vector<double> &&vec);

  vector(const vector &other);
  vector(vector &&other);

  double magnitude() const;
  double normalize();

  size_t size() const;

  double dot_product(const vector &other) const;

  vector operator+(const vector &other);
  vector operator-(const vector &other);

private:
  friend vector operator*(const vector &v, double f);
  friend vector operator*(double f, const vector &v);

  friend vector operator/(const vector &m, double f);

  std::vector<double> data_;
  size_t size;
};

} // namespace linalg