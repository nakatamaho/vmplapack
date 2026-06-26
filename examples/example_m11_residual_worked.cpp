// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"
#include "../tests/Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 90;
}

const char* status_name(vmplapack::Rstatus status) {
    switch (status) {
        case vmplapack::Rstatus::ok:
            return "ok";
        case vmplapack::Rstatus::unbounded:
            return "unbounded";
        case vmplapack::Rstatus::non_finite:
            return "non_finite";
        case vmplapack::Rstatus::invalid_input:
            return "invalid_input";
    }
    return "invalid_input";
}

mpfrxx::mpfr_class interval_midpoint(const vmplapack::oracle::Rdot_interval& interval,
                                     mpfr_prec_t precision) {
    mpfrxx::mpfr_class lo = vmplapack::oracle::widen_mpfr(interval.lo, precision);
    mpfrxx::mpfr_class hi = vmplapack::oracle::widen_mpfr(interval.hi, precision);
    mpfrxx::mpfr_class sum = mpfrxx::mpfr_class::with_precision(precision);
    sum = lo + hi;
    mpfrxx::mpfr_class half = mpfrxx::mpfr_class::with_precision(precision, 0.5);
    mpfrxx::mpfr_class mid = mpfrxx::mpfr_class::with_precision(precision);
    mid = sum * half;
    return mid;
}

template <class REAL>
bool display_covers_expected(const vmplapack::Rmidrad<REAL>& box, const REAL& expected) {
    REAL lo = box.mid - box.rad;
    REAL hi = box.mid + box.rad;
    return lo <= expected && expected <= hi;
}

template <class REAL>
void print_row(const std::vector<REAL>& matrix,
               const std::vector<REAL>& x,
               const std::vector<REAL>& b,
               const std::vector<vmplapack::Rmidrad<REAL>>& out,
               const std::vector<REAL>& expected,
               std::size_t row) {
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(x.size());
    const REAL* row_data = matrix.data() + row * static_cast<std::size_t>(n);
    vmplapack::oracle::Rdot_interval dot_ref = vmplapack::oracle::Rdot_oracle(n, row_data, 1, x.data(), 1);
    mpfr_prec_t precision = dot_ref.precision + static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits() + 128L);
    mpfrxx::mpfr_class dot_mid = interval_midpoint(dot_ref, precision);
    mpfrxx::mpfr_class b_mp = vmplapack::oracle::widen_value(b[row], precision);
    mpfrxx::mpfr_class residual_mid = b_mp - dot_mid;

    REAL lo = out[row].mid - out[row].rad;
    REAL hi = out[row].mid + out[row].rad;
    std::cout << "  row " << row << '\n';
    std::cout << "    oracle A[row]*x midpoint = " << dot_mid << '\n';
    std::cout << "    oracle residual midpoint = " << residual_mid << '\n';
    std::cout << "    expected residual = " << expected[row] << '\n';
    std::cout << "    mid = " << out[row].mid << '\n';
    std::cout << "    rad = " << out[row].rad << '\n';
    std::cout << "    lower display = " << lo << '\n';
    std::cout << "    upper display = " << hi << '\n';
    std::cout << "    status = " << status_name(out[row].status) << '\n';
    std::cout << "    expected covered by display = " << std::boolalpha
              << display_covers_expected(out[row], expected[row]) << '\n';
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL one = A::one();
    REAL two = one + one;
    REAL large = vmplapack::gendot::power_of_two<REAL>(exponent);
    REAL half = A::half();

    std::vector<REAL> x;
    x.push_back(one);
    x.push_back(one);
    x.push_back(one);

    std::vector<REAL> matrix;
    matrix.push_back(large);
    matrix.push_back(one);
    matrix.push_back(-large);
    matrix.push_back(one);
    matrix.push_back(-two);
    matrix.push_back(one);

    std::vector<REAL> b;
    b.push_back(one + half);
    b.push_back(-one);

    std::vector<REAL> expected;
    expected.push_back(half);
    expected.push_back(-one);

    std::vector<vmplapack::Rmidrad<REAL>> out(static_cast<std::size_t>(2));
    vmplapack::Rstatus status = vmplapack::vRresidual<REAL>(2, 3, matrix.data(), 3, x.data(), b.data(), out.data());

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  A = [[2^" << exponent << ", 1, -2^" << exponent << "], [1, -2, 1]]" << '\n';
    std::cout << "  x = [1, 1, 1]" << '\n';
    std::cout << "  b = [1.5, -1]" << '\n';
    std::cout << "  residual = b - A*x" << '\n';
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  return status_rank = " << vmplapack::Rstatus_rank(status) << '\n';
    print_row(matrix, x, b, out, expected, 0U);
    print_row(matrix, x, b, out, expected, 1U);
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return EXIT_SUCCESS;
}
