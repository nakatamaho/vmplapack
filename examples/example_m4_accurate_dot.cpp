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
    return 90;
}

template <class REAL>
REAL naive_dot(const REAL* x, const REAL* y, std::ptrdiff_t n) {
    using A = vmplapack::Rarith<REAL>;

    REAL acc = A::zero();
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        REAL prod = x[i] * y[i];
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
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

mpfrxx::mpfr_class upward_abs_diff(const mpfrxx::mpfr_class& a,
                                   const mpfrxx::mpfr_class& b,
                                   mpfr_prec_t precision) {
    mpfrxx::mpfr_class aa = vmplapack::oracle::widen_mpfr(a, precision);
    mpfrxx::mpfr_class bb = vmplapack::oracle::widen_mpfr(b, precision);
    mpfrxx::mpfr_class result = mpfrxx::mpfr_class::with_precision(precision);
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    if (aa >= bb) {
        result = aa - bb;
    } else {
        result = bb - aa;
    }
    return result;
}

mpfrxx::mpfr_class upward_condition_rhs(const vmplapack::oracle::Rdot_interval& interval,
                                        const mpfrxx::mpfr_class& ref,
                                        const mpfrxx::mpfr_class& s_hi,
                                        long tier_precision,
                                        mpfr_prec_t precision) {
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    mpfrxx::mpfr_class u = vmplapack::oracle::power_of_two(precision, -static_cast<mpfr_exp_t>(tier_precision));
    mpfrxx::mpfr_class abs_ref = vmplapack::oracle::abs_at(ref, precision);
    mpfrxx::mpfr_class term1 = mpfrxx::mpfr_class::with_precision(precision);
    term1 = u * abs_ref;

    mpfrxx::mpfr_class u2 = mpfrxx::mpfr_class::with_precision(precision);
    u2 = u * u;
    mpfrxx::mpfr_class c = mpfrxx::mpfr_class::with_precision(precision, 4.0);
    mpfrxx::mpfr_class cu2 = mpfrxx::mpfr_class::with_precision(precision);
    cu2 = c * u2;
    mpfrxx::mpfr_class term2 = mpfrxx::mpfr_class::with_precision(precision);
    term2 = cu2 * s_hi;

    mpfrxx::mpfr_class lo = vmplapack::oracle::widen_mpfr(interval.lo, precision);
    mpfrxx::mpfr_class hi = vmplapack::oracle::widen_mpfr(interval.hi, precision);
    mpfrxx::mpfr_class width = mpfrxx::mpfr_class::with_precision(precision);
    width = hi - lo;
    mpfrxx::mpfr_class floor_factor = mpfrxx::mpfr_class::with_precision(precision, 8.0);
    mpfrxx::mpfr_class floor = mpfrxx::mpfr_class::with_precision(precision);
    floor = floor_factor * width;

    mpfrxx::mpfr_class partial = mpfrxx::mpfr_class::with_precision(precision);
    partial = term1 + term2;
    mpfrxx::mpfr_class rhs = mpfrxx::mpfr_class::with_precision(precision);
    rhs = partial + floor;
    return rhs;
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    vmplapack::gendot::Rdot_case<REAL> c =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent, A::one());
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    REAL naive = naive_dot(c.x.data(), c.y.data(), n);
    REAL accurate = vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::oracle::Rdot_interval interval = vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    mpfr_prec_t precision = interval.precision + 128;
    mpfrxx::mpfr_class oracle_mid = interval_midpoint(interval, precision);
    mpfrxx::mpfr_class naive_mp = vmplapack::oracle::widen_value(naive, precision);
    mpfrxx::mpfr_class accurate_mp = vmplapack::oracle::widen_value(accurate, precision);
    mpfrxx::mpfr_class naive_error = upward_abs_diff(naive_mp, oracle_mid, precision);
    mpfrxx::mpfr_class rdot_error = upward_abs_diff(accurate_mp, oracle_mid, precision);
    mpfrxx::mpfr_class s_hi = vmplapack::oracle::upward_abs_term_sum(n, c.x.data(), 1, c.y.data(), 1, precision);
    mpfrxx::mpfr_class rhs = upward_condition_rhs(interval, oracle_mid, s_hi, A::precision_bits(), precision);
    mpfrxx::mpfr_class condition = vmplapack::oracle::condition_estimate(n,
                                                                         c.x.data(),
                                                                         1,
                                                                         c.y.data(),
                                                                         1,
                                                                         c.exact,
                                                                         precision);
    bool passes = rdot_error <= rhs;

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  case = [2^" << exponent << ", 1, 2^" << exponent << "] dot [1, 1, -1]" << '\n';
    std::cout << "  exact = " << c.exact << '\n';
    std::cout << "  naive dot = " << naive << '\n';
    std::cout << "  Rdot = " << accurate << '\n';
    std::cout << "  oracle midpoint = " << oracle_mid << '\n';
    std::cout << "  oracle lo = " << interval.lo << '\n';
    std::cout << "  oracle hi = " << interval.hi << '\n';
    std::cout << "  condition estimate = " << condition << '\n';
    std::cout << "  naive error = " << naive_error << '\n';
    std::cout << "  Rdot error = " << rdot_error << '\n';
    std::cout << "  M4 bound RHS = " << rhs << '\n';
    std::cout << "  Rdot passes bound = " << passes << '\n';
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return 0;
}
