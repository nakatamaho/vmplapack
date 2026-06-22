// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#ifndef VMPLAPACK_ENABLE_MPFR
#error "Rdot_oracle.h requires VMPLAPACK_ENABLE_MPFR. Native-only builds skip oracle tests."
#endif

#include <vmplapack/vmplapack.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace vmplapack {
namespace oracle {

struct Rdot_interval {
    mpfrxx::mpfr_class lo;
    mpfrxx::mpfr_class hi;
    mpfr_prec_t precision;
    int refinements;
};

inline mpfrxx::mpfr_class zero_at(mpfr_prec_t precision) {
    mpfrxx::mpfr_class z = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_set_zero(z.mpfr_data(), 0);
    return z;
}

inline mpfrxx::mpfr_class one_at(mpfr_prec_t precision) {
    mpfrxx::mpfr_class one = mpfrxx::mpfr_class::with_precision(precision, 1.0);
    return one;
}

inline mpfrxx::mpfr_class power_of_two(mpfr_prec_t precision, mpfr_exp_t exponent) {
    mpfrxx::mpfr_class value = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_set_ui_2exp(value.mpfr_data(), 1, exponent, MPFR_RNDN);
    return value;
}

template <class REAL>
mpfrxx::mpfr_class widen_value(const REAL& value, mpfr_prec_t precision) {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        mpfrxx::mpfr_class widened = mpfrxx::mpfr_class::with_precision(precision, static_cast<double>(value));
        return widened;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        mpfrxx::mpfr_class widened(value.mpfr_data(), precision);
        return widened;
    } else {
        static_assert(Ralways_false<REAL>::value, "oracle::widen_value does not support this REAL type.");
    }
}

template <class REAL>
bool widening_is_exact(const REAL& value, mpfr_prec_t precision) {
    mpfrxx::mpfr_class widened = widen_value(value, precision);
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return mpfr_cmp_d(widened.mpfr_data(), static_cast<double>(value)) == 0;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        return mpfr_get_prec(widened.mpfr_data()) == precision && mpfr_cmp(widened.mpfr_data(), value.mpfr_data()) == 0;
    } else {
        static_assert(Ralways_false<REAL>::value, "oracle::widening_is_exact does not support this REAL type.");
    }
}

template <class REAL>
mpfr_prec_t initial_precision() {
    long tier_precision = Rarith<REAL>::precision_bits();
    long heuristic = 4L * tier_precision + 64L;
    long chosen = std::max(512L, heuristic);
    return static_cast<mpfr_prec_t>(chosen);
}

inline mpfrxx::mpfr_class widen_mpfr(const mpfrxx::mpfr_class& value, mpfr_prec_t precision) {
    mpfrxx::mpfr_class widened(value.mpfr_data(), precision);
    return widened;
}

inline mpfrxx::mpfr_class abs_at(const mpfrxx::mpfr_class& value, mpfr_prec_t precision) {
    mpfrxx::mpfr_class widened = widen_mpfr(value, precision);
    mpfrxx::mpfr_class result = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_abs(result.mpfr_data(), widened.mpfr_data(), MPFR_RNDN);
    return result;
}

inline mpfrxx::mpfr_class min_at(const mpfrxx::mpfr_class& a, const mpfrxx::mpfr_class& b, mpfr_prec_t precision) {
    mpfrxx::mpfr_class aa = widen_mpfr(a, precision);
    mpfrxx::mpfr_class bb = widen_mpfr(b, precision);
    return (aa <= bb) ? aa : bb;
}

inline mpfrxx::mpfr_class max_at(const mpfrxx::mpfr_class& a, const mpfrxx::mpfr_class& b, mpfr_prec_t precision) {
    mpfrxx::mpfr_class aa = widen_mpfr(a, precision);
    mpfrxx::mpfr_class bb = widen_mpfr(b, precision);
    return (aa >= bb) ? aa : bb;
}

inline bool contains(const Rdot_interval& interval, const mpfrxx::mpfr_class& value) {
    mpfr_prec_t precision = std::max(interval.precision, mpfr_get_prec(value.mpfr_data()));
    mpfrxx::mpfr_class lo = widen_mpfr(interval.lo, precision);
    mpfrxx::mpfr_class hi = widen_mpfr(interval.hi, precision);
    mpfrxx::mpfr_class exact = widen_mpfr(value, precision);
    return lo <= exact && exact <= hi;
}

template <class REAL>
bool contains_exact(const Rdot_interval& interval, const REAL& value) {
    mpfr_prec_t precision = std::max(interval.precision, static_cast<mpfr_prec_t>(Rarith<REAL>::precision_bits() + 64L));
    mpfrxx::mpfr_class exact = widen_value(value, precision);
    return contains(interval, exact);
}

inline bool intervals_overlap(const Rdot_interval& a, const Rdot_interval& b) {
    mpfr_prec_t precision = std::max(a.precision, b.precision);
    mpfrxx::mpfr_class alo = widen_mpfr(a.lo, precision);
    mpfrxx::mpfr_class ahi = widen_mpfr(a.hi, precision);
    mpfrxx::mpfr_class blo = widen_mpfr(b.lo, precision);
    mpfrxx::mpfr_class bhi = widen_mpfr(b.hi, precision);
    return alo <= bhi && blo <= ahi;
}

inline bool intervals_agree_for_tier(const Rdot_interval& a, const Rdot_interval& b, long tier_precision) {
    if (!intervals_overlap(a, b)) {
        return false;
    }

    mpfr_prec_t precision = std::max(a.precision, b.precision) + 64;
    mpfrxx::mpfr_class lo_min = min_at(a.lo, b.lo, precision);
    mpfrxx::mpfr_class hi_max = max_at(a.hi, b.hi, precision);
    mpfrxx::mpfr_class width = hi_max - lo_min;

    mpfrxx::mpfr_class one = one_at(precision);
    mpfrxx::mpfr_class scale = one;
    mpfrxx::mpfr_class abs_lo = abs_at(lo_min, precision);
    mpfrxx::mpfr_class abs_hi = abs_at(hi_max, precision);
    if (abs_lo > scale) {
        scale = abs_lo;
    }
    if (abs_hi > scale) {
        scale = abs_hi;
    }

    mpfrxx::mpfr_class tolerance_unit = power_of_two(precision, -static_cast<mpfr_exp_t>(tier_precision + 16L));
    mpfrxx::mpfr_class tolerance = scale * tolerance_unit;
    return width <= tolerance;
}

template <class REAL>
Rdot_interval dot_interval_at_precision(std::ptrdiff_t n,
                                         const REAL* x,
                                         std::ptrdiff_t incx,
                                         const REAL* y,
                                         std::ptrdiff_t incy,
                                         mpfr_prec_t precision) {
    mpfrxx::mpfr_class lo = zero_at(precision);
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDD);
        mpfrxx::mpfr_class acc = zero_at(precision);
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            mpfrxx::mpfr_class xi = widen_value(x[i * incx], precision);
            mpfrxx::mpfr_class yi = widen_value(y[i * incy], precision);
            mpfrxx::mpfr_class prod = mpfrxx::mpfr_class::with_precision(precision);
            prod = xi * yi;
            mpfrxx::mpfr_class next = mpfrxx::mpfr_class::with_precision(precision);
            next = acc + prod;
            acc = next;
        }
        lo = acc;
    }

    mpfrxx::mpfr_class hi = zero_at(precision);
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
        mpfrxx::mpfr_class acc = zero_at(precision);
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            mpfrxx::mpfr_class xi = widen_value(x[i * incx], precision);
            mpfrxx::mpfr_class yi = widen_value(y[i * incy], precision);
            mpfrxx::mpfr_class prod = mpfrxx::mpfr_class::with_precision(precision);
            prod = xi * yi;
            mpfrxx::mpfr_class next = mpfrxx::mpfr_class::with_precision(precision);
            next = acc + prod;
            acc = next;
        }
        hi = acc;
    }

    return {lo, hi, precision, 0};
}

template <class REAL>
Rdot_interval Rdot_oracle(std::ptrdiff_t n,
                          const REAL* x,
                          std::ptrdiff_t incx,
                          const REAL* y,
                          std::ptrdiff_t incy,
                          int max_refinements = 8) {
    mpfr_prec_t precision = initial_precision<REAL>();
    for (int refinement = 0; refinement <= max_refinements; ++refinement) {
        Rdot_interval current = dot_interval_at_precision(n, x, incx, y, incy, precision);
        Rdot_interval doubled = dot_interval_at_precision(n, x, incx, y, incy, precision * 2);
        if (intervals_agree_for_tier(current, doubled, Rarith<REAL>::precision_bits())) {
            doubled.refinements = refinement;
            return doubled;
        }
        precision *= 2;
    }
    throw std::runtime_error("Rdot_oracle: precision doubling did not stabilize");
}

template <class REAL>
mpfrxx::mpfr_class upward_abs_term_sum(std::ptrdiff_t n,
                                       const REAL* x,
                                       std::ptrdiff_t incx,
                                       const REAL* y,
                                       std::ptrdiff_t incy,
                                       mpfr_prec_t precision) {
    mpfrxx::mpfr_class acc = zero_at(precision);
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        mpfrxx::mpfr_class xi = widen_value(x[i * incx], precision);
        mpfrxx::mpfr_class yi = widen_value(y[i * incy], precision);
        mpfrxx::mpfr_class ax = abs_at(xi, precision);
        mpfrxx::mpfr_class ay = abs_at(yi, precision);
        mpfrxx::mpfr_class prod = mpfrxx::mpfr_class::with_precision(precision);
        prod = ax * ay;
        mpfrxx::mpfr_class next = mpfrxx::mpfr_class::with_precision(precision);
        next = acc + prod;
        acc = next;
    }
    return acc;
}

template <class REAL>
mpfrxx::mpfr_class condition_estimate(std::ptrdiff_t n,
                                      const REAL* x,
                                      std::ptrdiff_t incx,
                                      const REAL* y,
                                      std::ptrdiff_t incy,
                                      const REAL& exact,
                                      mpfr_prec_t precision) {
    mpfrxx::mpfr_class sum_abs = upward_abs_term_sum(n, x, incx, y, incy, precision);
    mpfrxx::mpfr_class exact_mp = widen_value(exact, precision);
    mpfrxx::mpfr_class abs_exact = abs_at(exact_mp, precision);
    mpfrxx::mpfr_class two = mpfrxx::mpfr_class::with_precision(precision, 2.0);
    mpfrxx::mpfr_class numerator = two * sum_abs;
    mpfrxx::mpfr_class condition = numerator / abs_exact;
    return condition;
}

} // namespace oracle
} // namespace vmplapack
