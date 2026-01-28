#pragma once
#include <cstddef>
#include <vector>
#include <stdexcept>

class Matrix
{
public:
    Matrix(size_t rows, size_t cols);

    Matrix(const std::vector<std::vector<double>> &values);

    Matrix(size_t rows, size_t cols, const std::vector<double> &values);

    double &operator()(size_t rows, size_t cols);
    double operator()(size_t rows, size_t cols) const;

    bool operator==(const Matrix &other) const;
    bool operator!=(const Matrix &other) const;

    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;

    size_t rows() const;
    size_t cols() const;

    void rref();

private:
    void swap_rows(size_t r1, size_t r2);
    void scale_row(size_t r, double scalar);
    void add_row_multiple(size_t src, size_t dst, double scalar);

    friend std::ostream &operator<<(std::ostream &os, const Matrix &m);

    friend Matrix operator*(const Matrix &m, double f);
    friend Matrix operator*(double f, const Matrix &m);

    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};