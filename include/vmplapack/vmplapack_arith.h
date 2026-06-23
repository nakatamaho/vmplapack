// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cfenv>
#include <cmath>
#include <limits>

#ifdef VMPLAPACK_ENABLE_MPFR
#include <mpfrxx_mkII.h>
#endif

namespace vmplapack {

template <class REAL>
struct Rarith;

template <>
struct Rarith<float> {
    static long precision_bits() noexcept { return 24L; }
    static float unit_roundoff() noexcept { return 0x1p-24f; }
    static float zero() noexcept { return 0.0f; }
    static float one() noexcept { return 1.0f; }
    static float half() noexcept { return 0.5f; }
    static float infinity() noexcept { return std::numeric_limits<float>::infinity(); }
    static float eta() noexcept { return std::numeric_limits<float>::denorm_min(); }
    static bool is_finite(float x) noexcept { return std::isfinite(x); }
    static float abs(float x) noexcept { return std::fabs(x); }
    static float fma(float a, float b, float c) noexcept { return std::fmaf(a, b, c); }

    struct round_up {
        int saved_;

        round_up() noexcept : saved_(std::fegetround()) { (void)std::fesetround(FE_UPWARD); }
        ~round_up() noexcept { (void)std::fesetround(saved_); }

        round_up(const round_up&) = delete;
        round_up& operator=(const round_up&) = delete;
        round_up(round_up&&) = delete;
        round_up& operator=(round_up&&) = delete;
    };

    struct round_down {
        int saved_;

        round_down() noexcept : saved_(std::fegetround()) { (void)std::fesetround(FE_DOWNWARD); }
        ~round_down() noexcept { (void)std::fesetround(saved_); }

        round_down(const round_down&) = delete;
        round_down& operator=(const round_down&) = delete;
        round_down(round_down&&) = delete;
        round_down& operator=(round_down&&) = delete;
    };
};

template <>
struct Rarith<double> {
    static long precision_bits() noexcept { return 53L; }
    static double unit_roundoff() noexcept { return 0x1p-53; }
    static double zero() noexcept { return 0.0; }
    static double one() noexcept { return 1.0; }
    static double half() noexcept { return 0.5; }
    static double infinity() noexcept { return std::numeric_limits<double>::infinity(); }
    static double eta() noexcept { return std::numeric_limits<double>::denorm_min(); }
    static bool is_finite(double x) noexcept { return std::isfinite(x); }
    static double abs(double x) noexcept { return std::fabs(x); }
    static double fma(double a, double b, double c) noexcept { return std::fma(a, b, c); }

    struct round_up {
        int saved_;

        round_up() noexcept : saved_(std::fegetround()) { (void)std::fesetround(FE_UPWARD); }
        ~round_up() noexcept { (void)std::fesetround(saved_); }

        round_up(const round_up&) = delete;
        round_up& operator=(const round_up&) = delete;
        round_up(round_up&&) = delete;
        round_up& operator=(round_up&&) = delete;
    };

    struct round_down {
        int saved_;

        round_down() noexcept : saved_(std::fegetround()) { (void)std::fesetround(FE_DOWNWARD); }
        ~round_down() noexcept { (void)std::fesetround(saved_); }

        round_down(const round_down&) = delete;
        round_down& operator=(const round_down&) = delete;
        round_down(round_down&&) = delete;
        round_down& operator=(round_down&&) = delete;
    };
};

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
struct Rarith<mpfrxx::mpfr_class> {
    static long precision_bits() noexcept { return static_cast<long>(mpfrxx::default_prec()); }

    static mpfrxx::mpfr_class unit_roundoff() {
        mpfr_prec_t w = static_cast<mpfr_prec_t>(precision_bits());
        mpfrxx::mpfr_class u = mpfrxx::mpfr_class::with_precision(w);
        mpfr_exp_t exponent = -static_cast<mpfr_exp_t>(precision_bits());
        mpfr_set_ui_2exp(u.mpfr_data(), 1, exponent, MPFR_RNDN);
        return u;
    }

    static mpfrxx::mpfr_class zero() {
        return mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(precision_bits()), 0.0);
    }

    static mpfrxx::mpfr_class one() {
        return mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(precision_bits()), 1.0);
    }

    static mpfrxx::mpfr_class half() {
        return mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(precision_bits()), 0.5);
    }

    static mpfrxx::mpfr_class infinity() {
        mpfrxx::mpfr_class inf = mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(precision_bits()));
        mpfr_set_inf(inf.mpfr_data(), 1);
        return inf;
    }

    static mpfrxx::mpfr_class eta() {
        mpfrxx::mpfr_class value =
            mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(precision_bits()));
        mpfr_exp_t exponent = mpfr_get_emin() - 1;
        mpfr_set_ui_2exp(value.mpfr_data(), 1, exponent, MPFR_RNDN);
        return value;
    }

    static bool is_finite(const mpfrxx::mpfr_class& x) noexcept { return mpfr_number_p(x.mpfr_data()) != 0; }

    static mpfrxx::mpfr_class abs(const mpfrxx::mpfr_class& x) {
        mpfrxx::mpfr_class y = mpfrxx::abs(x);
        return y;
    }

    static mpfrxx::mpfr_class fma(const mpfrxx::mpfr_class& a,
                                  const mpfrxx::mpfr_class& b,
                                  const mpfrxx::mpfr_class& c) {
        mpfrxx::mpfr_class y = mpfrxx::fma(a, b, c);
        return y;
    }

    struct round_up {
        mpfrxx::rounding_mode_scope scope_;

        round_up() noexcept : scope_(MPFR_RNDU) {}

        round_up(const round_up&) = delete;
        round_up& operator=(const round_up&) = delete;
        round_up(round_up&&) = delete;
        round_up& operator=(round_up&&) = delete;
    };

    struct round_down {
        mpfrxx::rounding_mode_scope scope_;

        round_down() noexcept : scope_(MPFR_RNDD) {}

        round_down(const round_down&) = delete;
        round_down& operator=(const round_down&) = delete;
        round_down(round_down&&) = delete;
        round_down& operator=(round_down&&) = delete;
    };
};
#endif

} // namespace vmplapack
