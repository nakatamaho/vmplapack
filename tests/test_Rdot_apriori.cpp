// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"
#include "Rgendot.h"

#include <vmplapack/vmplapack.h>

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

template <class REAL>
REAL quiet_nan_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::quiet_NaN();
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_set_nan(value.mpfr_data());
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "quiet_nan_value does not support this REAL type.");
    }
}

template <class REAL>
REAL overflow_large_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::max();
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_exp_t exponent = mpfr_get_emax() - 1;
        mpfr_set_ui_2exp(value.mpfr_data(), 1, exponent, MPFR_RNDN);
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "overflow_large_value does not support this REAL type.");
    }
}

template <class REAL>
bool midrad_covers_oracle(const vmplapack::Rmidrad<REAL>& m,
                          const vmplapack::oracle::Rdot_interval& ref) {
    using A = vmplapack::Rarith<REAL>;

    if (m.status == vmplapack::Rstatus::unbounded) {
        return !A::is_finite(m.rad);
    }
    if (m.status != vmplapack::Rstatus::ok || m.rad < A::zero()) {
        return false;
    }

    mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 128L);
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(m.mid, precision);
    mpfrxx::mpfr_class rad = vmplapack::oracle::widen_value(m.rad, precision);
    mpfrxx::mpfr_class lo = mpfrxx::mpfr_class::with_precision(precision);
    mpfrxx::mpfr_class hi = mpfrxx::mpfr_class::with_precision(precision);
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDD);
        lo = mid - rad;
    }
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
        hi = mid + rad;
    }
    mpfrxx::mpfr_class ref_lo = vmplapack::oracle::widen_mpfr(ref.lo, precision);
    mpfrxx::mpfr_class ref_hi = vmplapack::oracle::widen_mpfr(ref.hi, precision);
    return lo <= ref_lo && ref_hi <= hi;
}

template <class REAL>
void require_abs_sum_upper_bounds_oracle(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    // Catches an S_up implementation that calls accurate Rdot or otherwise underestimates sum |x_i*y_i|.
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    REAL s_up = vmplapack::Rsum_abs_dot_upper(n, c.x.data(), 1, c.y.data(), 1);
    mpfr_prec_t precision = static_cast<mpfr_prec_t>(4L * vmplapack::Rarith<REAL>::precision_bits() + 512L);
    mpfrxx::mpfr_class s_up_mp = vmplapack::oracle::widen_value(s_up, precision);
    mpfrxx::mpfr_class oracle_s_up =
        vmplapack::oracle::upward_abs_term_sum(n, c.x.data(), 1, c.y.data(), 1, precision);
    require<REAL>(s_up_mp >= oracle_s_up, tier, "Rsum_abs_dot_upper underestimated the oracle upper sum");
}

template <class REAL>
void require_apriori_contains_oracle(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    // Catches a non-outward or wrong-constant a-priori radius by checking oracle interval inclusion.
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    vmplapack::Rmidrad<REAL> box = vmplapack::vRdot_apriori(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::oracle::Rdot_interval ref =
        vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    require<REAL>(box.status == vmplapack::Rstatus::ok, tier, "vRdot_apriori returned non-ok for a regular case");
    require<REAL>(box.mid == vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1),
                  tier,
                  "vRdot_apriori midpoint is not Dot2");
    require<REAL>(midrad_covers_oracle(box, ref), tier, "vRdot_apriori did not cover the oracle interval");
}

template <class REAL>
void test_regular_cases(const char* tier, const int* exponents, std::size_t exponent_count) {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL eps = (A::one() + (two * A::unit_roundoff())) - A::one();
    REAL y2 = -A::one() + eps;

    vmplapack::gendot::Rdot_case<REAL> family_a = vmplapack::gendot::family_a_alternating<REAL>(8, eps);
    require_abs_sum_upper_bounds_oracle(family_a, tier);
    require_apriori_contains_oracle(family_a, tier);

    vmplapack::gendot::Rdot_case<REAL> family_b = vmplapack::gendot::family_b_two_term<REAL>(y2);
    require_abs_sum_upper_bounds_oracle(family_b, tier);
    require_apriori_contains_oracle(family_b, tier);

    for (std::size_t i = 0; i < exponent_count; ++i) {
        vmplapack::gendot::Rdot_case<REAL> family_c =
            vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponents[i], A::one());
        require_abs_sum_upper_bounds_oracle(family_c, tier);
        require_apriori_contains_oracle(family_c, tier);
    }
}

template <class REAL>
void test_apriori_tighter_than_directed(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    // Catches the M5 directed interval being returned from the M7 API instead of the a-priori Dot2 bound.
    vmplapack::gendot::Rdot_case<REAL> c =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent, A::one());
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    vmplapack::Rmidrad<REAL> directed = vmplapack::vRdot(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::Rmidrad<REAL> apriori = vmplapack::vRdot_apriori(n, c.x.data(), 1, c.y.data(), 1);
    require<REAL>(directed.status == vmplapack::Rstatus::ok, tier, "directed comparison case returned non-ok");
    require<REAL>(apriori.status == vmplapack::Rstatus::ok, tier, "a-priori comparison case returned non-ok");
    require<REAL>(apriori.rad < directed.rad, tier, "a-priori radius was not tighter than the directed radius");
}

template <class REAL>
void test_strides_and_boundaries(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches a-priori routines that validate strides but then index as if both strides were one.
    std::vector<REAL> xs(static_cast<std::size_t>(10), A::zero());
    std::vector<REAL> ys(static_cast<std::size_t>(13), A::zero());
    xs[0] = A::one();
    xs[2] = -A::one();
    xs[4] = A::half();
    xs[6] = A::one() + A::one();
    ys[0] = A::one();
    ys[3] = A::one();
    ys[6] = -A::one();
    ys[9] = A::half();

    vmplapack::Rmidrad<REAL> strided = vmplapack::vRdot_apriori(4, xs.data(), 2, ys.data(), 3);
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(4, xs.data(), 2, ys.data(), 3);
    require<REAL>(strided.status == vmplapack::Rstatus::ok, tier, "strided vRdot_apriori returned non-ok");
    require<REAL>(midrad_covers_oracle(strided, ref), tier, "strided vRdot_apriori did not cover oracle interval");

    // Catches verified empty/invalid input behavior that accidentally follows the accurate UB contract.
    vmplapack::Rmidrad<REAL> empty = vmplapack::vRdot_apriori<REAL>(0, nullptr, 1, nullptr, 1);
    require<REAL>(empty.status == vmplapack::Rstatus::ok && empty.mid == A::zero() && empty.rad == A::zero(),
                  tier,
                  "vRdot_apriori n==0 did not return the exact zero enclosure");
    REAL one = A::one();
    require<REAL>(vmplapack::vRdot_apriori<REAL>(-1, nullptr, 1, nullptr, 1).status ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRdot_apriori n<0 was not invalid_input");
    require<REAL>(vmplapack::vRdot_apriori<REAL>(1, nullptr, 1, &one, 1).status ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRdot_apriori null x pointer was not invalid_input");
    require<REAL>(vmplapack::vRdot_apriori<REAL>(1, &one, 0, &one, 1).status ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRdot_apriori non-positive stride was not invalid_input");
}

template <class REAL>
void test_nonfinite_overflow_and_underflow(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches false finite certificates for NaN/Inf inputs.
    REAL nan_value = quiet_nan_value<REAL>();
    REAL one = A::one();
    require<REAL>(vmplapack::vRdot_apriori<REAL>(1, &nan_value, 1, &one, 1).status ==
                      vmplapack::Rstatus::non_finite,
                  tier,
                  "vRdot_apriori non-finite input was not non_finite");

    // Catches finite-input overflow reported as ok instead of unbounded.
    REAL large = overflow_large_value<REAL>();
    REAL two = A::one() + A::one();
    vmplapack::Rmidrad<REAL> overflow = vmplapack::vRdot_apriori<REAL>(1, &large, 1, &two, 1);
    require<REAL>(overflow.status == vmplapack::Rstatus::unbounded, tier, "vRdot_apriori overflow was not unbounded");
    require<REAL>(!A::is_finite(overflow.rad), tier, "vRdot_apriori unbounded overflow did not have infinite radius");

    // Catches a native underflow term that is missing or rounded inward.
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        REAL eta = A::eta();
        REAL half = A::half();
        vmplapack::Rmidrad<REAL> underflow = vmplapack::vRdot_apriori<REAL>(1, &eta, 1, &half, 1);
        vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(1, &eta, 1, &half, 1);
        require<REAL>(underflow.status == vmplapack::Rstatus::ok, tier, "native underflow case was not ok");
        require<REAL>(midrad_covers_oracle(underflow, ref), tier, "native underflow case did not cover oracle interval");
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL eta = A::eta();
        REAL half = A::half();
        vmplapack::Rmidrad<REAL> underflow = vmplapack::vRdot_apriori<REAL>(1, &eta, 1, &half, 1);
        require<REAL>(underflow.status == vmplapack::Rstatus::unbounded,
                      tier,
                      "MPFR underflow flag did not force an unbounded result");
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "unsupported REAL type");
    }
}

void test_integer_gate() {
    // Catches evaluating n*u<1 in floating point or shifting by the signed length type width.
    require<float>(vmplapack::detail::apriori_length_gate(7, 3), "gate", "n=7,p=3 should pass");
    require<float>(!vmplapack::detail::apriori_length_gate(8, 3), "gate", "n=8,p=3 should fail");
    require<float>(!vmplapack::detail::apriori_length_gate(-1, 24), "gate", "negative n should fail");
    long digits = static_cast<long>(std::numeric_limits<std::ptrdiff_t>::digits);
    require<float>(vmplapack::detail::apriori_length_gate(std::numeric_limits<std::ptrdiff_t>::max(), digits),
                   "gate",
                   "p==digits should pass without forming 2^p");
}

template <class REAL>
void test_tier(const char* tier, const int* exponents, std::size_t exponent_count) {
    test_regular_cases<REAL>(tier, exponents, exponent_count);
    test_apriori_tighter_than_directed<REAL>(tier, exponents[exponent_count - 1]);
    test_strides_and_boundaries<REAL>(tier);
    test_nonfinite_overflow_and_underflow<REAL>(tier);
}

} // namespace

int main() {
    long expected_mpfr_precision = parse_expected_mpfr_precision();
    if (expected_mpfr_precision != 0) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_mpfr_precision));
    }

    test_integer_gate();

    const int float_exponents[] = {24, 30};
    const int double_exponents[] = {53, 60};
    test_tier<float>("float", float_exponents, 2);
    test_tier<double>("double", double_exponents, 2);

    if (expected_mpfr_precision == 53) {
        test_tier<mpfrxx::mpfr_class>("mpfr", double_exponents, 2);
    } else if (expected_mpfr_precision == 512) {
        const int mpfr_w512_exponents[] = {520, 600, 700};
        test_tier<mpfrxx::mpfr_class>("mpfr", mpfr_w512_exponents, 3);
    } else if (expected_mpfr_precision != 0) {
        std::cerr << "unsupported MPFRXX_DEFAULT_PRECISION_BITS=" << expected_mpfr_precision << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
