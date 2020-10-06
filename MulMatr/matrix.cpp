//
// Created by alex on 27.09.2020.
//

#include <cassert>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <cmath>
#include "matrix.h"


matrix::matrix(size_t n) : n(n) {
    assert(n > 0);

    buf = new double[n * n];
}

matrix::matrix(matrix &m) : n(m.n) {
    buf = new double[n * n];
    memcpy(buf, m.buf, n * n * sizeof(double));
}

matrix::matrix(matrix &&m) noexcept : n(m.n), buf(m.buf) {
    m.buf = nullptr;
}

matrix &matrix::operator=(matrix const &m) {
    if (&m == this) return *this;

    n = m.n;
    double *tmp = new double[n * n];
    memcpy(tmp, m.buf, n * n * sizeof(double));
    delete[]buf;
    buf = tmp;

    return *this;
}

matrix &matrix::operator=(matrix &&m) noexcept {
    if (&m == this) return *this;

    n = m.n;
    delete[]buf;
    buf = m.buf;
    m.buf = nullptr;

    return *this;
}

matrix::~matrix() {
    delete[]buf;
}

matrix::proxy_row matrix::operator[](const size_t row) const {
    assert(row < n);

    const proxy_row pr(&buf[row * n], n);
    return pr;
}

matrix::proxy_row matrix::operator[](size_t row) {
    assert(row < n);

    proxy_row pr(&buf[row * n], n);
    return pr;
}

matrix matrix::operator+(const matrix &m) const {
    assert(n == m.n);
    matrix res(n);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            res[i][j] = (*this)[i][j] + m[i][j];
        }
    }

    return res;
}

matrix matrix::operator-(const matrix &m) const {
    assert(n == m.n);
    matrix res(n);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            res[i][j] = (*this)[i][j] - m[i][j];
        }
    }

    return res;
}

matrix matrix::operator*(const matrix &m) const {
    assert(n == m.n);
    matrix res(n);
    bzero(res.buf, n * n * sizeof(double));

    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                res[i][j] += (*this)[i][k] * m[k][j];
            }
        }
    }

    return res;
}

std::vector<double> matrix::operator*(const std::vector<double> &v) const {
    assert(v.size() == n);
    std::vector<double> res(n, 0.0);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            res[i] += (*this)[i][j] * v[j];
        }
    }

    return res;
}

double matrix::proxy_row::operator[](const size_t col) const {
    assert(col < n);

    return buf[col];
}

double &matrix::proxy_row::operator[](const size_t col) {
    assert(col < n);

    return buf[col];
}

std::ostream &operator<<(std::ostream &f, const matrix &m) {
    for (size_t i = 0; i < m.size(); i++) {
        for (size_t j = 0; j < m.size(); j++) {
            f << m[i][j] << ' ';
        }
        f << '\n';
    }
    f << '\n';
    return f;
}

size_t matrix::size() const {
    return n;
}

matrix matrix::transpose() const {
    matrix res(n);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            res[j][i] = (*this)[i][j];
        }
    }

    return res;
}

matrix matrix::get_block(size_t X, size_t Y) const {
    assert((X >= 1) && (X <= 2) && (Y >= 1) && (Y <= 2));
    size_t N_2 = n >> 1;
    matrix res(N_2);

    size_t i0 = (X - 1) * N_2;
    size_t j0 = (Y - 1) * N_2;
    for (size_t i = 0; i < N_2; i++) {
        for (size_t j = 0; j < N_2; j++) {
            res[i][j] = (*this)[i0 + i][j0 + j];
        }
    }
    return res;
}

void matrix::set_block(const matrix &m, size_t X, size_t Y) {
    assert((X >= 1) && (X <= 2) && (Y >= 1) && (Y <= 2));
    size_t N_2 = m.n;
    assert(n == (N_2 << 1));
    size_t i0 = (X - 1) * N_2;
    size_t j0 = (Y - 1) * N_2;
    for (size_t i = 0; i < N_2; i++) {
        for (size_t j = 0; j < N_2; j++) {
            (*this)[i0 + i][j0 + j] = m[i][j];
        }
    }
}

matrix::matrix(const matrix &m11, const matrix &m12, const matrix &m21, const matrix &m22) {
    assert((m11.n == m12.n) && (m12.n == m21.n) && (m21.n == m22.n));

    size_t N_2 = m11.n;

    n = N_2 << 1;
    buf = new double[n * n];

    set_block(m11, 1, 1);
    set_block(m12, 1, 2);
    set_block(m21, 2, 1);
    set_block(m22, 2, 2);
}

bool matrix::operator==(const matrix &m) const {
    if (n != m.n) return false;
    if (this == &m) return true;

    const double eps = 1e-8;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (std::abs((*this)[i][j] - m[i][j]) > eps) return false;
        }
    }

    return true;
}

matrix shtrassen_mul(const matrix &A, const matrix &B) {
    assert(A.size() == B.size());

    size_t N = A.size();
    size_t N_2 = N >> 1;

    if ((N < 128) || ((N & 1) == 1)) {
        return A * B;
    }

    matrix A11 = A.get_block(1, 1), A12 = A.get_block(1, 2), A21 = A.get_block(2, 1), A22 = A.get_block(2, 2);
    matrix B11 = B.get_block(1, 1), B12 = B.get_block(1, 2), B21 = B.get_block(2, 1), B22 = B.get_block(2, 2);

    matrix P1 = shtrassen_mul(A11 + A22, B11 + B22);
    matrix P2 = shtrassen_mul(A21 + A22, B11);
    matrix P3 = shtrassen_mul(A11, B12 - B22);
    matrix P4 = shtrassen_mul(A22, B21 - B11);
    matrix P5 = shtrassen_mul(A11 + A12, B22);
    matrix P6 = shtrassen_mul(A21 - A11, B11 + B12);
    matrix P7 = shtrassen_mul(A12 - A22, B21 + B22);

    matrix C11 = P1 + P4 - P5 + P7;
    matrix C12 = P3 + P5;
    matrix C21 = P2 + P4;
    matrix C22 = P1 - P2 + P3 + P6;

    matrix res(C11, C12, C21, C22);

    return res;
}
