//
// Created by alex on 27.09.2020.
//

#ifndef MULMATR_MATRIX_H
#define MULMATR_MATRIX_H


#include <cstddef>
#include <iostream>
#include <vector>

class matrix {
    class proxy_row;
public:
    matrix() = default;
    explicit matrix(size_t n);
    matrix(matrix& m);
    matrix(matrix&& m) noexcept;
    matrix& operator=(matrix const & m);
    matrix& operator=(matrix&& m) noexcept;
    ~matrix();

    matrix(const matrix& m11, const matrix& m12, const matrix& m21, const matrix& m22);

    proxy_row operator[](size_t row) const;
    proxy_row operator[](size_t row);

    matrix operator+(const matrix& m) const;
    matrix operator-(const matrix& m) const ;
    matrix operator*(const matrix& m) const;

    std::vector<double> operator*(const std::vector<double>& v) const;

    bool operator==(const matrix& m) const;

    matrix get_block(size_t X, size_t Y) const;

    void set_block(const matrix& m, size_t X, size_t Y);

    matrix transpose() const;

    size_t size() const;
private:
    size_t n = 0;
    double* buf=nullptr;

    class proxy_row {
    public:
        proxy_row(double* buf, size_t n) : n(n),buf(buf) {};
        double operator[](size_t col) const;
        double& operator[](size_t col);
    private:
        double* buf;
        size_t n;
    };
};

std::ostream& operator<<(std::ostream& fs, const matrix& m);

matrix shtrassen_mul(const matrix& A, const matrix& B);

#endif //MULMATR_MATRIX_H
