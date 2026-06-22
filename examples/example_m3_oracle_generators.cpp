// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"
#include "../tests/Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 180;
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
REAL exact_epsilon_above_one() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    REAL eps = next - A::one();
    return eps;
}

template <class REAL>
void show_terms(const vmplapack::gendot::Rdot_case<REAL>& c) {
    using A = vmplapack::Rarith<REAL>;

    REAL prefix = A::zero();
    std::cout << "  terms" << '\n';
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        REAL product = c.x[i] * c.y[i];
        REAL next = prefix + product;
        prefix = next;
        std::cout << "    i=" << i
                  << " x=" << c.x[i]
                  << " y=" << c.y[i]
                  << " x*y=" << product
                  << " naive_prefix=" << prefix << '\n';
    }
}

template <class REAL>
void show_case(const char* tier,
               const char* name,
               const vmplapack::gendot::Rdot_case<REAL>& c) {
    mpfr_prec_t precision = vmplapack::oracle::initial_precision<REAL>() * 2;
    vmplapack::oracle::Rdot_interval interval =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(c.x.size()), c.x.data(), 1, c.y.data(), 1);
    mpfrxx::mpfr_class condition =
        vmplapack::oracle::condition_estimate(static_cast<std::ptrdiff_t>(c.x.size()),
                                              c.x.data(),
                                              1,
                                              c.y.data(),
                                              1,
                                              c.exact,
                                              precision);
    mpfrxx::mpfr_class mid = interval_midpoint(interval);

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << tier << " / " << name << '\n';
    std::cout << "  n = " << c.x.size() << '\n';
    show_terms(c);
    std::cout << "  expected exact = " << c.exact << '\n';
    std::cout << "  oracle lo = " << interval.lo << '\n';
    std::cout << "  oracle hi = " << interval.hi << '\n';
    std::cout << "  oracle midpoint = " << mid << '\n';
    std::cout << "  oracle precision = " << interval.precision << '\n';
    std::cout << "  refinements = " << interval.refinements << '\n';
    std::cout << "  condition estimate = " << condition << '\n';
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL eps = exact_epsilon_above_one<REAL>();
    REAL y2 = -A::one() + eps;
    vmplapack::gendot::Rdot_case<REAL> family_a = vmplapack::gendot::family_a_alternating<REAL>(4, eps);
    vmplapack::gendot::Rdot_case<REAL> family_b = vmplapack::gendot::family_b_two_term<REAL>(y2);
    vmplapack::gendot::Rdot_case<REAL> family_c =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent, A::one());

    show_case(tier, "Family A", family_a);
    show_case(tier, "Family B", family_b);
    show_case(tier, "Family C", family_c);
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return 0;
}
