// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"
#include "../tests/Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstddef>
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

template <class REAL>
REAL naive_sum(const std::vector<REAL>& x) {
    using A = vmplapack::Rarith<REAL>;

    REAL acc = A::zero();
    for (std::size_t i = 0; i < x.size(); ++i) {
        REAL next = acc + x[i];
        acc = next;
    }
    return acc;
}

template <class REAL>
std::vector<REAL> ones(std::size_t n) {
    std::vector<REAL> y;
    y.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        y.push_back(vmplapack::Rarith<REAL>::one());
    }
    return y;
}

mpfrxx::mpfr_class interval_midpoint(const vmplapack::oracle::Rdot_interval& interval) {
    mpfr_prec_t precision = interval.precision + 64;
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
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL large = vmplapack::gendot::power_of_two<REAL>(exponent);
    std::vector<REAL> x;
    x.push_back(large);
    x.push_back(A::one());
    x.push_back(-large);
    std::vector<REAL> y = ones<REAL>(x.size());

    REAL naive = naive_sum(x);
    REAL accurate = vmplapack::Rsum(static_cast<std::ptrdiff_t>(x.size()), x.data(), 1);
    vmplapack::oracle::Rdot_interval interval =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(x.size()), x.data(), 1, y.data(), 1);
    mpfrxx::mpfr_class oracle_mid = interval_midpoint(interval);

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << tier << '\n';
    std::cout << "  case = [2^" << exponent << ", 1, -2^" << exponent << "]" << '\n';
    std::cout << "  exact = " << A::one() << '\n';
    std::cout << "  naive sum = " << naive << '\n';
    std::cout << "  Rsum = " << accurate << '\n';
    std::cout << "  oracle midpoint = " << oracle_mid << '\n';
    std::cout << "  oracle lo = " << interval.lo << '\n';
    std::cout << "  oracle hi = " << interval.hi << '\n';
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return 0;
}
