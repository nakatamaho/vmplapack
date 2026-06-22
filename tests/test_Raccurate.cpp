// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"
#include "Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

long parse_expected_mpfr_precision() {
    const char* text = std::getenv("MPFRXX_DEFAULT_PRECISION_BITS");
    if (text == nullptr || text[0] == '\0') {
        return 0;
    }
    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0) {
        std::cerr << "invalid MPFRXX_DEFAULT_PRECISION_BITS=" << text << '\n';
        std::exit(EXIT_FAILURE);
    }
    return value;
}

mpfrxx::mpfr_class midpoint(const vmplapack::oracle::Rdot_interval& interval, mpfr_prec_t precision) {
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
REAL exact_epsilon_above_one() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    REAL eps = next - A::one();
    return eps;
}

template <class REAL>
std::vector<REAL> ones(std::size_t n) {
    std::vector<REAL> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        result.push_back(vmplapack::Rarith<REAL>::one());
    }
    return result;
}

template <class REAL>
void require_mpfr_case_precision(const vmplapack::gendot::Rdot_case<REAL>&, long) {}

template <>
void require_mpfr_case_precision<mpfrxx::mpfr_class>(const vmplapack::gendot::Rdot_case<mpfrxx::mpfr_class>& c,
                                                     long precision) {
    require<mpfrxx::mpfr_class>(vmplapack::gendot::all_inputs_have_precision(c, precision),
                                "mpfr",
                                "generator inputs are not at W");
}

template <class REAL>
void require_rdot_condition_bound(std::ptrdiff_t n,
                                  const REAL* x,
                                  std::ptrdiff_t incx,
                                  const REAL* y,
                                  std::ptrdiff_t incy,
                                  const char* tier) {
    // Catches a non-compensated dot that is only judged by a relative-to-reference tolerance.
    REAL got = vmplapack::Rdot(n, x, incx, y, incy);
    vmplapack::oracle::Rdot_interval interval = vmplapack::oracle::Rdot_oracle(n, x, incx, y, incy);
    mpfr_prec_t precision = interval.precision + 128;
    mpfrxx::mpfr_class got_mp = vmplapack::oracle::widen_value(got, precision);
    mpfrxx::mpfr_class ref = midpoint(interval, precision);
    mpfrxx::mpfr_class diff = upward_abs_diff(got_mp, ref, precision);
    mpfrxx::mpfr_class s_hi = vmplapack::oracle::upward_abs_term_sum(n, x, incx, y, incy, precision);
    mpfrxx::mpfr_class rhs = upward_condition_rhs(interval,
                                                  ref,
                                                  s_hi,
                                                  vmplapack::Rarith<REAL>::precision_bits(),
                                                  precision);
    require<REAL>(diff <= rhs, tier, "Rdot failed the condition-aware Dot2 bound");
}

template <class REAL>
void require_rdot_condition_bound(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    require_rdot_condition_bound(static_cast<std::ptrdiff_t>(c.x.size()), c.x.data(), 1, c.y.data(), 1, tier);
}

template <class REAL>
void require_rsum_condition_bound(std::ptrdiff_t n,
                                  const REAL* x,
                                  std::ptrdiff_t incx,
                                  const char* tier) {
    // Catches a Sum2 regression by testing sum as dot(x, 1) against the same interval oracle.
    std::vector<REAL> y = ones<REAL>(static_cast<std::size_t>(n));
    REAL got = vmplapack::Rsum(n, x, incx);
    vmplapack::oracle::Rdot_interval interval = vmplapack::oracle::Rdot_oracle(n, x, incx, y.data(), 1);
    mpfr_prec_t precision = interval.precision + 128;
    mpfrxx::mpfr_class got_mp = vmplapack::oracle::widen_value(got, precision);
    mpfrxx::mpfr_class ref = midpoint(interval, precision);
    mpfrxx::mpfr_class diff = upward_abs_diff(got_mp, ref, precision);
    mpfrxx::mpfr_class s_hi = vmplapack::oracle::upward_abs_term_sum(n, x, incx, y.data(), 1, precision);
    mpfrxx::mpfr_class rhs = upward_condition_rhs(interval,
                                                  ref,
                                                  s_hi,
                                                  vmplapack::Rarith<REAL>::precision_bits(),
                                                  precision);
    require<REAL>(diff <= rhs, tier, "Rsum failed the condition-aware Sum2 bound");
}

template <class REAL>
void require_rsum_condition_bound(const std::vector<REAL>& x, const char* tier) {
    require_rsum_condition_bound(static_cast<std::ptrdiff_t>(x.size()), x.data(), 1, tier);
}

template <class REAL>
std::vector<REAL> alternating_sum_case(int pairs, REAL delta) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> x;
    x.reserve(static_cast<std::size_t>(2 * pairs + 1));
    for (int i = 0; i < pairs; ++i) {
        x.push_back(A::one());
        x.push_back(-A::one());
    }
    x.push_back(delta);
    return x;
}

template <class REAL>
std::vector<REAL> exponent_sum_case(int exponent) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> x;
    REAL large = vmplapack::gendot::power_of_two<REAL>(exponent);
    x.push_back(large);
    x.push_back(A::one());
    x.push_back(-large);
    return x;
}

template <class REAL>
void test_accurate_boundaries(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches accurate routines that inspect null pointers for empty reductions.
    require<REAL>(vmplapack::Rsum<REAL>(0, nullptr, 1) == A::zero(), tier, "Rsum n==0 did not return zero");
    require<REAL>(vmplapack::Rsum<REAL>(-3, nullptr, 7) == A::zero(), tier, "Rsum n<0 did not return zero");
    require<REAL>(vmplapack::Rdot<REAL>(0, nullptr, 1, nullptr, 1) == A::zero(), tier, "Rdot n==0 did not return zero");
    require<REAL>(vmplapack::Rdot<REAL>(-2, nullptr, 7, nullptr, 9) == A::zero(), tier, "Rdot n<0 did not return zero");
}

template <class REAL>
void test_strided_cases(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches incorrect BLAS stride indexing in accurate sum and dot.
    std::vector<REAL> xs(static_cast<std::size_t>(13), A::zero());
    std::vector<REAL> ys(static_cast<std::size_t>(10), A::zero());
    xs[0] = A::one();
    xs[3] = -A::one();
    xs[6] = exact_epsilon_above_one<REAL>();
    xs[9] = A::half();
    ys[0] = A::one();
    ys[2] = A::one();
    ys[4] = A::one();
    ys[6] = -A::one();
    require_rsum_condition_bound(static_cast<std::ptrdiff_t>(4), xs.data(), 3, tier);
    require_rdot_condition_bound(static_cast<std::ptrdiff_t>(4), xs.data(), 3, ys.data(), 2, tier);
}

template <class REAL>
void test_accurate_tier(const char* tier, const int* exponents, std::size_t exponent_count, long mpfr_precision) {
    using A = vmplapack::Rarith<REAL>;

    test_accurate_boundaries<REAL>(tier);

    REAL eps = exact_epsilon_above_one<REAL>();
    REAL y2 = -A::one() + eps;

    vmplapack::gendot::Rdot_case<REAL> family_a = vmplapack::gendot::family_a_alternating<REAL>(8, eps);
    vmplapack::gendot::Rdot_case<REAL> family_b = vmplapack::gendot::family_b_two_term<REAL>(y2);
    require_mpfr_case_precision(family_a, mpfr_precision);
    require_mpfr_case_precision(family_b, mpfr_precision);
    require_rdot_condition_bound(family_a, tier);
    require_rdot_condition_bound(family_b, tier);
    require_rsum_condition_bound(alternating_sum_case<REAL>(8, eps), tier);
    require_rsum_condition_bound(std::vector<REAL>{A::one(), y2}, tier);

    for (std::size_t i = 0; i < exponent_count; ++i) {
        vmplapack::gendot::Rdot_case<REAL> c =
            vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponents[i], A::one());
        require_mpfr_case_precision(c, mpfr_precision);
        require_rdot_condition_bound(c, tier);
        require_rsum_condition_bound(exponent_sum_case<REAL>(exponents[i]), tier);
    }

    test_strided_cases<REAL>(tier);
}

void test_cross_tier_double_mpfr_w53() {
    using REAL = mpfrxx::mpfr_class;

    // Catches hidden double detours or MPFR materialization errors under the binary64-like W=53 tier.
    const double xd[] = {1.25, -2.0, 0.5, 8.0, -0.125};
    const double yd[] = {3.0, -0.25, 2.0, -1.0, 6.0};
    std::vector<REAL> xm;
    std::vector<REAL> ym;
    for (std::size_t i = 0; i < 5; ++i) {
        xm.push_back(REAL::with_precision(53, xd[i]));
        ym.push_back(REAL::with_precision(53, yd[i]));
    }

    double got_double = vmplapack::Rdot<double>(5, xd, 1, yd, 1);
    REAL got_mpfr = vmplapack::Rdot<REAL>(5, xm.data(), 1, ym.data(), 1);

    mpfr_prec_t precision = 256;
    mpfrxx::mpfr_class d = vmplapack::oracle::widen_value(got_double, precision);
    mpfrxx::mpfr_class m = vmplapack::oracle::widen_value(got_mpfr, precision);
    mpfrxx::mpfr_class diff = upward_abs_diff(d, m, precision);
    double mag = std::fabs(got_double);
    double ulp = std::nextafter(mag, std::numeric_limits<double>::infinity()) - mag;
    mpfrxx::mpfr_class tolerance = mpfrxx::mpfr_class::with_precision(precision, 2.0 * ulp);
    require<REAL>(diff <= tolerance, "mpfr", "Rdot<double> and Rdot<mpfr@53> differed by more than two result ulps");
}

} // namespace

int main() {
    long expected_mpfr_precision = parse_expected_mpfr_precision();
    if (expected_mpfr_precision != 0) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_mpfr_precision));
    }

    const int float_exponents[] = {24, 30};
    const int double_exponents[] = {53, 60};
    test_accurate_tier<float>("float", float_exponents, 2, 0);
    test_accurate_tier<double>("double", double_exponents, 2, 0);

    if (expected_mpfr_precision == 53) {
        test_accurate_tier<mpfrxx::mpfr_class>("mpfr", double_exponents, 2, expected_mpfr_precision);
        test_cross_tier_double_mpfr_w53();
    } else if (expected_mpfr_precision == 512) {
        const int mpfr_w512_exponents[] = {520, 600, 700};
        test_accurate_tier<mpfrxx::mpfr_class>("mpfr", mpfr_w512_exponents, 3, expected_mpfr_precision);
    } else if (expected_mpfr_precision != 0) {
        std::cerr << "unsupported MPFRXX_DEFAULT_PRECISION_BITS=" << expected_mpfr_precision << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
